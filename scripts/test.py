#!/usr/bin/env pvpython

import argparse
import sys

import paraview
from paraview.simple import (  # type: ignore
    ColorBy,
    GetActiveViewOrCreate,
    GetColorTransferFunction,
    GetOpacityTransferFunction,
    GetScalarBar,
    Hide,
    Render,
    SaveScreenshot,
    Show,
    Threshold,
    XMLUnstructuredGridReader,
    _DisableFirstRenderCameraReset,
)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Visualize DisPerSE walls (.vtu) in ParaView/pvpython."
    )
    parser.add_argument(
        "--input",
        default="outputs/snap_010/snap_010_walls.vtu",
        help="Path to the VTU file exported by analyze_snapshot.py.",
    )
    parser.add_argument(
        "--field",
        default="field_value",
        help="Point-data field used for coloring (e.g., field_value, mass, log_field_value).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        nargs=2,
        metavar=("MIN", "MAX"),
        help="Optional lower/upper bounds in the chosen field; outside values are discarded.",
    )
    parser.add_argument(
        "--invert-colormap",
        action="store_true",
        help="Flip the active color transfer function.",
    )
    parser.add_argument(
        "--screenshot",
        help="If provided, save a PNG at this path instead of keeping the interactive window.",
    )
    parser.add_argument(
        "--resolution",
        type=int,
        nargs=2,
        default=(1920, 1080),
        metavar=("W", "H"),
        help="Screenshot resolution in pixels.",
    )
    parser.add_argument(
        "--scalar-bar-position",
        type=float,
        nargs=2,
        default=(0.87, 0.12),
        metavar=("X", "Y"),
        help="Lower-left corner of the scalar bar in normalized view coords.",
    )
    parser.add_argument(
        "--scalar-bar-size",
        type=float,
        nargs=2,
        default=(0.05, 0.25),
        metavar=("W", "H"),
        help="Scalar bar width/height in normalized view coords.",
    )
    parser.add_argument(
        "--scalar-bar-font-sizes",
        type=int,
        nargs=2,
        default=(12, 10),
        metavar=("TITLE", "LABEL"),
        help="Font sizes for the scalar bar title and tick labels.",
    )
    return parser.parse_args(argv)


def _get_array_info(data_info, array_name: str):
    """Compatibility helper to extract vtkPVArrayInformation by name."""
    getter = getattr(data_info, "GetArrayInformation", None)
    array_info = None
    if callable(getter):
        try:
            array_info = getter(array_name)
        except TypeError:
            array_info = None
        if array_info is None:
            num = getattr(data_info, "GetNumberOfArrays", lambda: 0)()
            for idx in range(num):
                candidate = getter(idx)
                name = getattr(candidate, "GetName", lambda: None)()
                if name == array_name:
                    array_info = candidate
                    break
    if array_info is None:
        getter = getattr(data_info, "GetArray", None)
        if callable(getter):
            try:
                array_info = getter(array_name)
            except TypeError:
                array_info = None
    return array_info


def _get_data_range(source, array_name: str, association: str) -> tuple[float, float]:
    """Return the min/max of `array_name` on the provided ParaView source."""
    info = source.GetDataInformation()
    if association == "POINTS":
        data_info = info.GetPointDataInformation()
    else:
        data_info = info.GetCellDataInformation()
    array_info = _get_array_info(data_info, array_name)
    if array_info is None:
        raise RuntimeError(
            f"Array '{array_name}' not found on {association.lower()} of the dataset."
        )
    get_range = getattr(array_info, "GetRange", None)
    if callable(get_range):
        rng = get_range()
    else:
        rng = getattr(array_info, "GetComponentRange", lambda *_: None)(-1)
    if rng is None:
        raise RuntimeError(f"Unable to determine range for array '{array_name}'.")
    return rng


def build_pipeline(args: argparse.Namespace):
    _DisableFirstRenderCameraReset()

    reader = XMLUnstructuredGridReader(FileName=[args.input])
    reader.PointArrayStatus = [args.field]

    source = reader
    if args.threshold is not None:
        thresh = Threshold(Input=reader)
        thresh.Scalars = ["POINTS", args.field]
        thresh.LowerThreshold, thresh.UpperThreshold = args.threshold
        source = thresh

    view = GetActiveViewOrCreate("RenderView")
    display = Show(source, view)
    display.SetRepresentationType("Surface")
    ColorBy(display, ("POINTS", args.field))
    display.SetScalarBarVisibility(view, True)

    # Older ParaView builds (e.g., 6.0.x) lack the VectorMode property that newer
    # RescaleTransferFunctionToDataRange() helpers touch, so avoid them and rescale
    # manually using the server-side data information instead.
    source.UpdatePipeline()
    rng = _get_data_range(source, args.field, "POINTS")

    lut = GetColorTransferFunction(args.field)
    pwf = GetOpacityTransferFunction(args.field)
    if lut is not None:
        lut.RescaleTransferFunction(rng[0], rng[1])
        if args.invert_colormap:
            lut.InvertTransferFunction()
        scalar_bar = GetScalarBar(lut, view)
        size = [float(args.scalar_bar_size[0]), float(args.scalar_bar_size[1])]
        try:
            scalar_bar.WindowLocation = "Any Location"
        except (AttributeError, paraview.NotSupportedException):
            pass
        try:
            scalar_bar.Position = [
                float(args.scalar_bar_position[0]),
                float(args.scalar_bar_position[1]),
            ]
        except (AttributeError, paraview.NotSupportedException):
            pass
        try:
            scalar_bar.Position2 = size
        except (AttributeError, paraview.NotSupportedException):
            # ParaView <=5.4 exposes ScalarBarLength/Thickness instead of Position2.
            length = size[1]
            width = size[0]
            try:
                scalar_bar.ScalarBarLength = length
            except (AttributeError, paraview.NotSupportedException):
                pass
            try:
                thickness_px = max(8, int(round(width * 200)))
                scalar_bar.ScalarBarThickness = thickness_px
            except (AttributeError, paraview.NotSupportedException):
                pass
        try:
            scalar_bar.TitleFontSize = int(args.scalar_bar_font_sizes[0])
        except (AttributeError, paraview.NotSupportedException):
            pass
        try:
            scalar_bar.LabelFontSize = int(args.scalar_bar_font_sizes[1])
        except (AttributeError, paraview.NotSupportedException):
            pass
        for prop in ("DrawTickMarks", "DrawTickLabels"):
            try:
                setattr(scalar_bar, prop, 1)
            except (AttributeError, paraview.NotSupportedException):
                pass
    if pwf is not None:
        pwf.RescaleTransferFunction(rng[0], rng[1])

    # Hide the reader when showing a thresholded source to reduce duplicated geometry.
    if source is not reader:
        Hide(reader, view)

    view.ResetCamera()
    return view


def main(argv: list[str]) -> None:
    args = parse_args(argv)
    view = build_pipeline(args)

    if args.screenshot:
        SaveScreenshot(
            args.screenshot,
            view,
            ImageResolution=tuple(args.resolution),
            TransparentBackground=False,
        )
    else:
        Render()


if __name__ == "__main__":
    main(sys.argv[1:])
