import os
from datetime import datetime, timedelta
from textwrap import shorten
import logging

# can os determine the terminal size?
try:
    n_cols = os.get_terminal_size().columns
    DYNAMIC_TERMINAL_WIDTH = True
    LINEBREAK_CHAR = '\r'
except Exception:
    DYNAMIC_TERMINAL_WIDTH = False
    LINEBREAK_CHAR = '\n'

# set function get_terminal_width depending on DYNAMIC_TERMINAL_WIDTH
if DYNAMIC_TERMINAL_WIDTH:
    def get_terminal_width(min_: int, default: int) -> int:
        return max(min_, os.get_terminal_size().columns)
else:
    def get_terminal_width(min_: int, default: int) -> int:
        return default


def stringify_timedelta(delta: timedelta) -> str:
    """
    Returns string 5 characters long.
    """
    d = delta.days
    h, rem = divmod(delta.seconds, 3600)
    m, s = divmod(rem, 60)
    if d:
        res = f'{d}d {h}h'
    elif h:
        res = f'{h}h {m}m'
    elif m:
        res = f'{m}m {s}s'
    else:
        res = f'{s}s'
    res = shorten(res, width=5, placeholder='')
    return f'{res:>5s}' if res else f'>999d'


def print_progress(
        i: int,
        n: int,
        message: str,
        start_time: datetime,
        message_width: int = 40,
        default_width: int = 100,
        sep: str = ' | ',
        end: str = LINEBREAK_CHAR
) -> None:
    message = f"{shorten(message, width=message_width, placeholder='...'):{message_width}s}"
    assert len(message) == message_width

    n = max(1, n)
    i_safe = min(max(1, i), n)
    time_left = stringify_timedelta((datetime.now() - start_time) / i_safe * (n - i_safe))  # 5 chars
    percentage = f"{f'{i / n:.0%}':>4s}"  # 4 chars

    width_total = get_terminal_width(min_=message_width + len(sep) * 2 + 20, default=default_width)

    text = f'{percentage}{sep}{time_left}{sep}{message}'
    len_progressbar = width_total - len(text)
    n_bars = len_progressbar - 3  # because of '[] '

    res = f"[{'=' * round(i / n * n_bars):{n_bars}}] {text}"

    if not len(res) == width_total:
        logging.warning('Failed to print progressbar!')

    print(res, end=end)
