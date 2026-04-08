#!/usr/bin/env python3
"""Export Codex and Claude Code sessions to Markdown and HTML.

Scans:
  - sessions/codex/... (Codex rollout-*.jsonl files)
  - sessions/claude/   (Claude Code UUID-named .jsonl files)

Writes to sessions/export/:
  - One .md per session
  - INDEX.md
  - MASTER.md
  - MASTER_BY_DAY.md
  - MASTER_BY_DAY.html
"""
import html
import json
import os
import re
from collections import Counter
from datetime import datetime, timezone

ROOT = os.path.abspath(os.path.dirname(__file__))
EXPORT_DIR = os.path.join(ROOT, "export")
CLAUDE_DIR = os.path.join(ROOT, "claude")

CODEX_FILENAME_RE = re.compile(
    r"rollout-(\d{4}-\d{2}-\d{2})T(\d{2})-(\d{2})-(\d{2})-([0-9a-f-]+)\.jsonl$"
)
CLAUDE_FILENAME_RE = re.compile(
    r"^([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})\.jsonl$"
)

SYSTEM_TAG_RE = re.compile(
    r"<(?:ide_opened_file|system-reminder|ide_selection)[^>]*>.*?</(?:ide_opened_file|system-reminder|ide_selection)>",
    re.DOTALL,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_iso(ts):
    if not ts:
        return None
    try:
        dt = datetime.fromisoformat(ts.replace("Z", "+00:00"))
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return dt
    except Exception:
        return None


def iso_to_readable(ts):
    dt = parse_iso(ts)
    if not dt:
        return ts
    return dt.strftime("%Y-%m-%d %H:%M:%SZ")


def safe_text(text):
    if text is None:
        return ""
    if not isinstance(text, str):
        return str(text)
    return text.rstrip()


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def find_jsonl_files():
    """Find all session JSONL files: Codex rollout-*.jsonl and Claude UUID.jsonl."""
    jsonl_files = []
    for root, _, files in os.walk(ROOT):
        abs_root = os.path.abspath(root)
        if abs_root.startswith(EXPORT_DIR):
            continue
        for name in files:
            if not name.endswith(".jsonl"):
                continue
            if CODEX_FILENAME_RE.search(name) or CLAUDE_FILENAME_RE.match(name):
                jsonl_files.append(os.path.join(root, name))
    return sorted(jsonl_files)


# ---------------------------------------------------------------------------
# Codex format
# ---------------------------------------------------------------------------

def parse_codex_filename_info(path):
    filename = os.path.basename(path)
    match = CODEX_FILENAME_RE.search(filename)
    if match:
        date_str, hh, mm, ss, session_id = match.groups()
        time_str = f"{hh}:{mm}:{ss}"
        start_dt = datetime.fromisoformat(f"{date_str}T{hh}:{mm}:{ss}+00:00")
        return {
            "date": date_str,
            "time": time_str,
            "session_id": session_id,
            "start_dt": start_dt,
            "ai_tool": "codex",
        }
    return {
        "date": "unknown-date",
        "time": "unknown-time",
        "session_id": None,
        "start_dt": None,
        "ai_tool": "codex",
    }


def parse_codex_session(path):
    session_meta = None
    messages = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except Exception:
                continue
            if obj.get("type") == "session_meta":
                session_meta = obj.get("payload") or {}
                continue
            if obj.get("type") == "event_msg":
                payload = obj.get("payload") or {}
                msg_type = payload.get("type")
                if msg_type in ("user_message", "agent_message"):
                    role = "user" if msg_type == "user_message" else "assistant"
                    messages.append({
                        "timestamp": obj.get("timestamp"),
                        "role": role,
                        "text": safe_text(payload.get("message")),
                    })
    if session_meta is None:
        session_meta = {}
    session_meta["ai_tool"] = "codex"
    return session_meta, messages


# ---------------------------------------------------------------------------
# Claude Code format
# ---------------------------------------------------------------------------

def extract_claude_text(content):
    """Extract human-readable text from a Claude Code content array.

    - Concatenates all text blocks, stripping system context tags.
    - Tool-use blocks are noted as [tool: name].
    - Tool-result and other blocks are skipped.
    """
    if isinstance(content, str):
        text = SYSTEM_TAG_RE.sub("", content).strip()
        return text
    if not isinstance(content, list):
        return ""
    parts = []
    for block in content:
        if not isinstance(block, dict):
            continue
        btype = block.get("type", "")
        if btype == "text":
            text = block.get("text", "")
            text = SYSTEM_TAG_RE.sub("", text).strip()
            if text:
                parts.append(text)
        elif btype == "tool_use":
            name = block.get("name", "tool")
            parts.append(f"[{name}]")
        # tool_result, image, etc. — skip
    return "\n".join(parts)


def parse_claude_filename_info(path):
    """Derive date/time from the first timestamped message in the file."""
    filename = os.path.basename(path)
    match = CLAUDE_FILENAME_RE.match(filename)
    session_id = match.group(1) if match else os.path.splitext(filename)[0]

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except Exception:
                continue
            if obj.get("type") in ("user", "assistant"):
                ts = obj.get("timestamp")
                dt = parse_iso(ts)
                if dt:
                    return {
                        "date": dt.strftime("%Y-%m-%d"),
                        "time": dt.strftime("%H:%M:%S"),
                        "session_id": session_id,
                        "start_dt": dt,
                        "ai_tool": "claude-code",
                    }
    return {
        "date": "unknown-date",
        "time": "unknown-time",
        "session_id": session_id,
        "start_dt": None,
        "ai_tool": "claude-code",
    }


def parse_claude_session(path):
    """Parse a Claude Code JSONL file into (session_meta, messages)."""
    messages = []
    session_id = None
    cwd = None
    model = None

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except Exception:
                continue

            obj_type = obj.get("type")

            # Skip sidechain (tool orchestration) messages
            if obj.get("isSidechain"):
                continue

            if obj_type == "user":
                if not session_id:
                    session_id = obj.get("sessionId")
                if not cwd:
                    cwd = obj.get("cwd")
                content = obj.get("message", {}).get("content", [])
                text = extract_claude_text(content)
                if not text:
                    continue
                messages.append({
                    "timestamp": obj.get("timestamp"),
                    "role": "user",
                    "text": text,
                })

            elif obj_type == "assistant":
                content = obj.get("message", {}).get("content", [])
                text = extract_claude_text(content)
                if not model:
                    model = obj.get("message", {}).get("model")
                if not text:
                    continue
                messages.append({
                    "timestamp": obj.get("timestamp"),
                    "role": "assistant",
                    "text": text,
                })

    session_meta = {
        "id": session_id,
        "cwd": cwd,
        "model_provider": "anthropic",
        "ai_tool": "claude-code",
        "model": model,
    }
    return session_meta, messages


# ---------------------------------------------------------------------------
# Format detection and dispatch
# ---------------------------------------------------------------------------

def is_claude_code_file(path):
    return bool(CLAUDE_FILENAME_RE.match(os.path.basename(path)))


def parse_filename_info(path):
    if is_claude_code_file(path):
        return parse_claude_filename_info(path)
    return parse_codex_filename_info(path)


def parse_session(path):
    if is_claude_code_file(path):
        return parse_claude_session(path)
    return parse_codex_session(path)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def render_session(path, session_meta, messages, file_info):
    date_str = file_info["date"]
    time_str = file_info["time"]
    ai_tool = file_info.get("ai_tool", "unknown")
    session_id = file_info["session_id"] or (session_meta.get("id") if session_meta else "unknown-id")

    out_name = f"{date_str}__{time_str.replace(':','-')}__{session_id}.md"
    out_path = os.path.join(EXPORT_DIR, out_name)

    meta_lines = []
    meta_lines.append(f"- Source file: `{os.path.relpath(path, ROOT)}`")
    if session_id:
        meta_lines.append(f"- Session id: `{session_id}`")
    meta_lines.append(f"- AI tool: `{ai_tool}`")
    if session_meta:
        ts = session_meta.get("timestamp") or file_info.get("start_dt")
        if isinstance(ts, datetime):
            ts = ts.strftime("%Y-%m-%d %H:%M:%SZ")
        if ts:
            meta_lines.append(f"- Session timestamp: `{iso_to_readable(ts) if isinstance(ts, str) else ts}`")
        if session_meta.get("cwd"):
            meta_lines.append(f"- CWD: `{session_meta.get('cwd')}`")
        if session_meta.get("model"):
            meta_lines.append(f"- Model: `{session_meta.get('model')}`")
        git = session_meta.get("git") or {}
        if git.get("commit_hash"):
            meta_lines.append(f"- Git commit: `{git.get('commit_hash')}`")
        if git.get("branch"):
            meta_lines.append(f"- Git branch: `{git.get('branch')}`")
        if session_meta.get("model_provider"):
            meta_lines.append(f"- Model provider: `{session_meta.get('model_provider')}`")
        if session_meta.get("cli_version"):
            meta_lines.append(f"- CLI version: `{session_meta.get('cli_version')}`")
    meta_lines.append(f"- Messages: `{len(messages)}`")

    lines = []
    lines.append(f"# Session {date_str} {time_str}")
    lines.append("")
    lines.append("## Metadata")
    lines.extend(meta_lines)
    lines.append("")
    lines.append("## Conversation")

    for msg in messages:
        idx = msg.get("index") or 0
        ts = iso_to_readable(msg.get("timestamp")) or "unknown-time"
        role = msg.get("role") or "unknown"
        anchor = f"m-{idx:04d}"
        lines.append(f'<a id="{anchor}"></a>')
        lines.append(f"### {idx}. {role} — {ts}")
        lines.append("")
        text = msg.get("text") or ""
        lines.append(text if text else "_No content_")
        lines.append("")

    return out_path, "\n".join(lines) + "\n", {
        "date": date_str,
        "time": time_str,
        "session_id": session_id,
        "ai_tool": ai_tool,
        "out_name": out_name,
        "message_count": len(messages),
        "cwd": session_meta.get("cwd") if session_meta else None,
    }


def render_index(entries):
    lines = ["# AI Sessions Index", ""]
    current_date = None
    for entry in entries:
        if entry["date"] != current_date:
            current_date = entry["date"]
            lines.append(f"## {current_date}")
        tool = entry.get("ai_tool", "")
        tool_badge = f" [{tool}]" if tool else ""
        label = f"{entry['time']}{tool_badge} — `{entry['session_id']}` ({entry['message_count']} msgs)"
        lines.append(f"- [{label}]({entry['out_name']})")
    lines.append("")
    return "\n".join(lines)


def render_master(entries, all_messages):
    lines = ["# AI Sessions Master Timeline", ""]
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")
    lines.append(f"Generated: `{generated}`")
    lines.append("")
    lines.append("## Sessions")
    for entry in entries:
        tool = entry.get("ai_tool", "")
        tool_badge = f" [{tool}]" if tool else ""
        label = f"{entry['date']} {entry['time']}{tool_badge} — `{entry['session_id']}` ({entry['message_count']} msgs)"
        lines.append(f"- [{label}]({entry['out_name']})")
    lines.append("")
    lines.append("## Timeline")

    for idx, msg in enumerate(all_messages, start=1):
        anchor = f"t-{idx:06d}"
        dt = msg.get("ts_dt")
        ts_label = dt.strftime("%Y-%m-%d %H:%M:%SZ") if dt else "unknown-time"
        role = msg.get("role") or "unknown"
        session_id = msg.get("session_id") or "unknown-id"
        session_date = msg.get("session_date") or "unknown-date"
        session_time = msg.get("session_time") or "unknown-time"
        session_out = msg.get("session_out_name") or ""
        session_idx = msg.get("index") or 0
        session_link = f"{session_out}#m-{session_idx:04d}" if session_out else ""

        lines.append(f'<a id="{anchor}"></a>')
        lines.append(f"### {ts_label} — {role} — `{session_id}`")
        lines.append("")
        if session_link:
            lines.append(f"Session: [{session_date} {session_time} — `{session_id}`]({session_link})")
            lines.append("")
        text = msg.get("text") or ""
        lines.append(text if text else "_No content_")
        lines.append("")

    return "\n".join(lines) + "\n"


def excerpt(text, limit):
    text = safe_text(text).replace("\n", " ").strip()
    if len(text) <= limit:
        return text if text else "_No content_"
    return text[: limit - 1].rstrip() + "…"


def slugify(text):
    keep = []
    for ch in text.lower():
        if ch.isalnum():
            keep.append(ch)
        elif ch in ("-", "_"):
            keep.append(ch)
        else:
            keep.append("-")
    slug = re.sub(r"-{2,}", "-", "".join(keep)).strip("-")
    return slug or "item"


def format_time(dt, fallback="unknown-time"):
    if not dt:
        return fallback
    return dt.strftime("%H:%M:%SZ")


STOPWORDS = {
    "a", "about", "above", "after", "again", "against", "all", "am", "an", "and",
    "any", "are", "as", "at", "be", "because", "been", "before", "being", "below",
    "between", "both", "but", "by", "can", "did", "do", "does", "doing", "down",
    "during", "each", "few", "for", "from", "further", "had", "has", "have", "having",
    "he", "her", "here", "hers", "herself", "him", "himself", "his", "how", "i",
    "if", "in", "into", "is", "it", "its", "itself", "just", "me", "more", "most",
    "my", "myself", "no", "nor", "not", "now", "of", "off", "on", "once", "only",
    "or", "other", "our", "ours", "ourselves", "out", "over", "own", "s", "same",
    "she", "should", "so", "some", "such", "t", "than", "that", "the", "their",
    "theirs", "them", "themselves", "then", "there", "these", "they", "this", "those",
    "through", "to", "too", "under", "until", "up", "very", "was", "we", "were",
    "what", "when", "where", "which", "while", "who", "whom", "why", "will", "with",
    "you", "your", "yours", "yourself", "yourselves", "http", "https", "www", "com",
    "org", "json", "jsonl", "md", "html", "png", "jpg", "jpeg", "gif",
}


def shift_headers(text, shift=3):
    """Shift all markdown headers in text down by `shift` levels (max depth 6)."""
    def replacer(m):
        new_level = min(len(m.group(1)) + shift, 6)
        return "#" * new_level + m.group(2)
    return re.sub(r'^(#{1,6})([ \t])', replacer, text, flags=re.MULTILINE)


def extract_keywords(messages, limit=10):
    counts = Counter()
    for msg in messages:
        text = safe_text(msg.get("text")).lower()
        for token in re.findall(r"[a-z][a-z0-9_-]{2,}", text):
            if token in STOPWORDS or token.isdigit():
                continue
            counts[token] += 1
    return [word for word, _ in counts.most_common(limit)]


def render_master_by_day(all_messages):
    lines = ["# AI Log", ""]
    lines.append("A chronological record of all AI-assisted development sessions for this project, including interactions with Codex and Claude Code.")
    lines.append("")
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")
    lines.append(f"Generated: `{generated}`")
    lines.append("")

    day_sessions = {}
    day_order = []

    for msg in all_messages:
        dt = msg.get("ts_dt") or msg.get("session_start_dt")
        day_key = dt.strftime("%Y-%m-%d") if dt else "unknown-date"
        if day_key not in day_sessions:
            day_sessions[day_key] = {}
            day_order.append(day_key)
        session_key = msg.get("session_out_name") or msg.get("session_id") or "unknown-session"
        if session_key not in day_sessions[day_key]:
            day_sessions[day_key][session_key] = {
                "session_id": msg.get("session_id") or "unknown-id",
                "ai_tool": msg.get("ai_tool") or "",
                "session_out_name": msg.get("session_out_name") or "",
                "session_date": msg.get("session_date") or "unknown-date",
                "session_time": msg.get("session_time") or "unknown-time",
                "messages": [],
                "first_dt": msg.get("ts_dt") or msg.get("session_start_dt"),
            }
        session = day_sessions[day_key][session_key]
        if not session["first_dt"]:
            session["first_dt"] = msg.get("ts_dt") or msg.get("session_start_dt")
        session["messages"].append(msg)

    for day_key in day_order:
        lines.append(f"## {day_key}")
        lines.append("")
        sessions = list(day_sessions[day_key].values())
        sessions.sort(key=lambda s: s["first_dt"] or datetime.max.replace(tzinfo=timezone.utc))
        day_messages = []
        for s in sessions:
            day_messages.extend(s["messages"])
        total_msgs = len(day_messages)
        lines.append(f"Sessions: `{len(sessions)}`  Messages: `{total_msgs}`")
        lines.append("")
        lines.append("### Summary (keywords)")
        keywords = extract_keywords(day_messages, limit=12)
        lines.append("Top keywords: " + (", ".join(f"`{k}`" for k in keywords) if keywords else "_None_"))
        lines.append("")
        for session in sessions:
            session_id = session["session_id"]
            session_out = session["session_out_name"]
            tool = session.get("ai_tool", "")
            tool_badge = f" [{tool}]" if tool else ""
            time_label = format_time(session["first_dt"], session["session_time"])
            session_label = f"{time_label}{tool_badge} — `{session_id}` ({len(session['messages'])} msgs)"
            # Use Quarto callout div — Pandoc handles these correctly unlike <details>
            escaped_label = session_label.replace('"', "'")
            lines.append(f'::: {{.callout-note collapse="true" title="{escaped_label}"}}')
            lines.append("")
            for msg in session["messages"]:
                text = msg.get("text") or ""
                # skip messages that are only tool-call labels e.g. "[Grep]\n[Read]"
                if not re.sub(r'\[[A-Za-z_]+\]\s*', '', text).strip():
                    continue
                dt = msg.get("ts_dt") or msg.get("session_start_dt")
                time_label = dt.strftime("%H:%M:%SZ") if dt else "unknown-time"
                role = msg.get("role") or "unknown"
                lines.append(f"#### {time_label} — {role}")
                lines.append("")
                lines.append(shift_headers(text, shift=4))
                lines.append("")
            lines.append(":::")
            lines.append("")

    return "\n".join(lines) + "\n"


def render_master_by_day_html(all_messages):
    day_sessions = {}
    day_order = []

    for msg in all_messages:
        dt = msg.get("ts_dt") or msg.get("session_start_dt")
        day_key = dt.strftime("%Y-%m-%d") if dt else "unknown-date"
        if day_key not in day_sessions:
            day_sessions[day_key] = {}
            day_order.append(day_key)
        session_key = msg.get("session_out_name") or msg.get("session_id") or "unknown-session"
        if session_key not in day_sessions[day_key]:
            day_sessions[day_key][session_key] = {
                "session_id": msg.get("session_id") or "unknown-id",
                "ai_tool": msg.get("ai_tool") or "",
                "session_out_name": msg.get("session_out_name") or "",
                "session_date": msg.get("session_date") or "unknown-date",
                "session_time": msg.get("session_time") or "unknown-time",
                "messages": [],
                "first_dt": msg.get("ts_dt") or msg.get("session_start_dt"),
            }
        session = day_sessions[day_key][session_key]
        if not session["first_dt"]:
            session["first_dt"] = msg.get("ts_dt") or msg.get("session_start_dt")
        session["messages"].append(msg)

    sidebar_items = []
    main_sections = []
    for day_key in day_order:
        day_id = f"day-{slugify(day_key)}"
        sidebar_items.append((day_key, day_id, []))
        sessions = list(day_sessions[day_key].values())
        sessions.sort(key=lambda s: s["first_dt"] or datetime.max.replace(tzinfo=timezone.utc))
        for session in sessions:
            session_id = session["session_id"]
            session_anchor = f"session-{slugify(day_key)}-{slugify(session_id)}-{slugify(session['session_time'])}"
            sidebar_items[-1][2].append((session, session_anchor))

        day_messages = []
        for s in sessions:
            day_messages.extend(s["messages"])
        total_msgs = len(day_messages)
        keywords = extract_keywords(day_messages, limit=10)
        section_lines = []
        section_lines.append(f'<section id="{day_id}">')
        section_lines.append(f"<h2>{html.escape(day_key)}</h2>")
        section_lines.append(
            f'<div class="day-meta">Sessions: <strong>{len(sessions)}</strong>'
            f' &nbsp; Messages: <strong>{total_msgs}</strong></div>'
        )
        section_lines.append('<div class="day-summary"><h3>Summary (keywords)</h3>')
        if keywords:
            section_lines.append(
                '<div class="keywords">' + " · ".join(html.escape(k) for k in keywords) + "</div>"
            )
        else:
            section_lines.append('<div class="keywords">No keywords</div>')
        section_lines.append("</div>")

        for session in sessions:
            session_id = session["session_id"]
            session_out = session["session_out_name"]
            tool = session.get("ai_tool", "")
            tool_badge = f" [{tool}]" if tool else ""
            session_anchor = f"session-{slugify(day_key)}-{slugify(session_id)}-{slugify(session['session_time'])}"
            time_label = format_time(session["first_dt"], session["session_time"])
            session_label = f"{time_label}{tool_badge} — {session_id} ({len(session['messages'])} msgs)"
            section_lines.append(
                f'<details id="{session_anchor}" data-session-label="{html.escape(session_label)}">'
            )
            section_lines.append(f"<summary>{html.escape(session_label)}</summary>")
            if session_out:
                section_lines.append(
                    f'<div class="session-link">Session file: '
                    f'<a href="{html.escape(session_out)}">{html.escape(session_out)}</a></div>'
                )
            section_lines.append('<div class="messages">')
            for msg in session["messages"]:
                dt = msg.get("ts_dt") or msg.get("session_start_dt")
                time_label = dt.strftime("%H:%M:%SZ") if dt else "unknown-time"
                role = msg.get("role") or "unknown"
                idx = msg.get("index") or 0
                link = f'{session_out}#m-{idx:04d}' if session_out else ""
                text = msg.get("text") or ""
                text_html = html.escape(text).replace("\n", "<br>")
                data_text = html.escape(safe_text(text).lower())
                meta = f"{time_label} — {role}"
                if link:
                    meta = f'{meta} — <a href="{html.escape(link)}">link</a>'
                section_lines.append(
                    f'<div class="msg" data-text="{data_text}"><div class="meta">{meta}</div>'
                    f'<div class="text">{text_html or "<em>No content</em>"}</div></div>'
                )
            section_lines.append("</div>")
            section_lines.append("</details>")
        section_lines.append("</section>")
        main_sections.append("\n".join(section_lines))

    sidebar_html = ['<nav class="sidebar">', "<h2>Timeline</h2>"]
    for day_key, day_id, sessions in sidebar_items:
        sidebar_html.append(f'<div class="side-day"><a href="#{day_id}">{html.escape(day_key)}</a>')
        sidebar_html.append('<div class="side-sessions">')
        for session, session_anchor in sessions:
            tool = session.get("ai_tool", "")
            tool_badge = f" [{tool}]" if tool else ""
            time_label = format_time(session["first_dt"], session["session_time"])
            label = f"{time_label}{tool_badge} — {session['session_id']}"
            sidebar_html.append(f'<a href="#{session_anchor}">{html.escape(label)}</a>')
        sidebar_html.append("</div></div>")
    sidebar_html.append("</nav>")

    styles = """
    :root {
      --bg: #f4f5f7; --panel: #ffffff; --ink: #1b1f24; --muted: #6b7280;
      --accent: #0b5fff; --border: #d6d9de;
      --mono: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
      --sans: -apple-system, BlinkMacSystemFont, "Segoe UI", "Helvetica Neue", Arial, sans-serif;
    }
    body { margin: 0; background: var(--bg); color: var(--ink); font-family: var(--sans); }
    .layout { display: flex; gap: 24px; padding: 24px; }
    .sidebar {
      width: 280px; min-width: 240px; max-height: calc(100vh - 48px);
      position: sticky; top: 24px; align-self: flex-start; overflow: auto;
      background: var(--panel); border: 1px solid var(--border); padding: 16px;
    }
    .sidebar h2 { margin: 0 0 12px 0; font-size: 16px; }
    .side-day { margin-bottom: 12px; }
    .side-day > a { font-weight: 600; color: var(--accent); text-decoration: none; }
    .side-sessions { display: flex; flex-direction: column; gap: 6px; margin-top: 6px; }
    .side-sessions a { color: var(--muted); text-decoration: none; font-size: 13px; }
    main { flex: 1; display: flex; flex-direction: column; gap: 24px; }
    section { background: var(--panel); border: 1px solid var(--border); padding: 20px 24px; }
    h1, h2, h3 { margin: 0 0 12px 0; }
    h1 { font-size: 22px; } h2 { font-size: 18px; } h3 { font-size: 15px; }
    .day-meta { color: var(--muted); font-size: 14px; margin-bottom: 12px; }
    .day-summary { margin-bottom: 16px; }
    .keywords { font-size: 13px; color: var(--muted); }
    details { border: 1px solid var(--border); padding: 12px 14px; margin-bottom: 12px; background: #fbfbfc; }
    summary { cursor: pointer; font-weight: 600; }
    .session-link { margin: 8px 0 12px; font-size: 13px; color: var(--muted); }
    .messages { display: flex; flex-direction: column; gap: 12px; }
    .msg { border-left: 2px solid var(--accent); padding-left: 12px; }
    .meta { font-size: 12px; color: var(--muted); font-family: var(--mono); margin-bottom: 6px; }
    .text { white-space: pre-wrap; font-size: 15px; line-height: 1.45; }
    .search { display: flex; gap: 8px; margin-bottom: 12px; }
    .search input { flex: 1; padding: 8px 10px; border: 1px solid var(--border); background: #fff; font-size: 14px; }
    .search button { padding: 8px 10px; border: 1px solid var(--border); background: #fff; cursor: pointer; }
    @media (max-width: 980px) { .layout { flex-direction: column; } .sidebar { width: auto; position: static; max-height: none; } }
    """

    html_lines = [
        "<!doctype html>",
        '<html><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        "<title>AI Sessions Timeline</title>",
        f"<style>{styles}</style>",
        "</head><body>",
        '<div class="layout">',
        "\n".join(sidebar_html),
        "<main>",
        "<section>",
        "<h1>AI Sessions Timeline</h1>",
        '<div class="search"><input id="search" placeholder="Filter sessions or messages">',
        '<button id="clear">Clear</button></div>',
        f'<div class="day-meta">Generated: <span class="meta">'
        f'{html.escape(datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ"))}'
        f"</span></div>",
        "</section>",
        "\n".join(main_sections),
        "</main></div>",
        "<script>\n"
        "const search = document.getElementById('search');\n"
        "const clearBtn = document.getElementById('clear');\n"
        "function applyFilter() {\n"
        "  const q = (search.value || '').trim().toLowerCase();\n"
        "  document.querySelectorAll('details').forEach(d => {\n"
        "    const sessionLabel = (d.getAttribute('data-session-label') || '').toLowerCase();\n"
        "    const msgs = d.querySelectorAll('.msg');\n"
        "    let anyMatch = false;\n"
        "    const sessionMatch = q && sessionLabel.includes(q);\n"
        "    msgs.forEach(m => {\n"
        "      const match = q ? (m.getAttribute('data-text') || '').includes(q) : true;\n"
        "      m.style.display = (sessionMatch || match) ? '' : 'none';\n"
        "      if (sessionMatch || match) anyMatch = true;\n"
        "    });\n"
        "    d.style.display = (!q || anyMatch) ? '' : 'none';\n"
        "    if (q && anyMatch) d.open = true;\n"
        "    if (!q) d.open = false;\n"
        "  });\n"
        "}\n"
        "search.addEventListener('input', applyFilter);\n"
        "clearBtn.addEventListener('click', () => { search.value = ''; applyFilter(); });\n"
        "</script>\n"
        "</body></html>",
    ]
    return "\n".join(html_lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(EXPORT_DIR, exist_ok=True)

    # Archive stale unknown-date exports from previous runs
    bad_files = [
        name for name in os.listdir(EXPORT_DIR)
        if name.startswith("unknown-date__unknown-time__") and name.endswith(".md")
    ]
    if bad_files:
        archive_dir = os.path.join(EXPORT_DIR, "_previous_run")
        os.makedirs(archive_dir, exist_ok=True)
        for name in bad_files:
            os.replace(os.path.join(EXPORT_DIR, name), os.path.join(archive_dir, name))

    entries = []
    all_messages = []

    for path in find_jsonl_files():
        file_info = parse_filename_info(path)
        session_meta, messages = parse_session(path)
        if not messages:
            continue
        ai_tool = file_info.get("ai_tool", "unknown")
        for idx, msg in enumerate(messages, start=1):
            msg["index"] = idx
            msg["session_id"] = file_info["session_id"] or (
                session_meta.get("id") if session_meta else "unknown-id"
            )
            msg["ai_tool"] = ai_tool
            msg["session_date"] = file_info["date"]
            msg["session_time"] = file_info["time"]
            session_meta_dt = parse_iso(session_meta.get("timestamp")) if session_meta else None
            msg["ts_dt"] = (
                parse_iso(msg.get("timestamp"))
                or session_meta_dt
                or file_info["start_dt"]
            )
            msg["session_start_dt"] = file_info["start_dt"]

        out_path, content, entry = render_session(path, session_meta, messages, file_info)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(content)
        for msg in messages:
            msg["session_out_name"] = entry["out_name"]
        all_messages.extend(messages)
        entries.append(entry)

    entries.sort(key=lambda e: (e["date"], e["time"], e["session_id"]))

    with open(os.path.join(EXPORT_DIR, "INDEX.md"), "w", encoding="utf-8") as f:
        f.write(render_index(entries))

    def sort_key(msg):
        dt = msg.get("ts_dt") or msg.get("session_start_dt")
        return (dt or datetime.max.replace(tzinfo=timezone.utc), msg.get("session_id") or "", msg.get("index") or 0)

    all_messages.sort(key=sort_key)

    with open(os.path.join(EXPORT_DIR, "MASTER.md"), "w", encoding="utf-8") as f:
        f.write(render_master(entries, all_messages))

    with open(os.path.join(EXPORT_DIR, "MASTER_BY_DAY.md"), "w", encoding="utf-8") as f:
        f.write(render_master_by_day(all_messages))

    with open(os.path.join(EXPORT_DIR, "MASTER_BY_DAY.html"), "w", encoding="utf-8") as f:
        f.write(render_master_by_day_html(all_messages))

    print(f"Wrote {len(entries)} sessions to {os.path.relpath(EXPORT_DIR, ROOT)}/")
    codex_count = sum(1 for e in entries if e.get("ai_tool") == "codex")
    claude_count = sum(1 for e in entries if e.get("ai_tool") == "claude-code")
    print(f"  Codex: {codex_count}  Claude Code: {claude_count}")


if __name__ == "__main__":
    main()
