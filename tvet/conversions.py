import datetime

MONTHS = {
    "jan": 1, "feb": 2, "mar": 3, "apr": 4,
    "may": 5, "jun": 6, "jul": 7, "aug": 8,
    "sep": 9, "oct": 10, "nov": 11, "dec": 12
}

def iso_to_jd(iso: str) -> float:
    iso = iso.strip()

    if not iso:
        raise ValueError("Datetime string is empty.")

    try: 
        if "T" in iso:
            date, time = iso.split("T", 1)
        else:
            date, time = iso.split(None, 1)

        year_str, month_str, day_str = date.split("-")
        hour_str, minute_str = time.split(":", 1)

        year = int(year_str)
        day = int(day_str)
        hour = int(hour_str)
        minute = int(minute_str)

    except ValueError as exc:
        raise ValueError(
            "Invalid datetime format. Expected "
            "'YYYY-MM-DDTHH:MM', 'YYYY-Mon-DDTHH:MM', "
            "'YYYY-MM-DD HH:MM', or 'YYYY-Mon-DD HH:MM'."
        ) from exc

    if month_str.isdigit():
        month = int(month_str)
    else:
        key = month_str.strip().lower()[:3]

        if key not in MONTHS:
            raise ValueError(f"Invalid month: {month_str}.")
        
        month = MONTHS[key]

    if (year, month, day) < (1582, 10, 15):
        raise ValueError("Dates before 1582-10-15 are not supported.")

    try:
        dt = datetime.datetime(year, month, day, hour, minute)
    except ValueError as exc:
        raise ValueError(f"Invalid calendar date/time: {exc}") from exc

    a = (14 - dt.month) // 12
    y = dt.year + 4800 - a
    m = dt.month + 12 * a - 3

    jdn = (
        dt.day
        + ((153 * m + 2) // 5)
        + 365 * y
        + y // 4
        - y // 100
        + y // 400
        - 32045
    )

    jd = jdn + (dt.hour - 12) / 24 + dt.minute / 1440

    return jd

def jd_to_iso(jd: float) -> str:
    if not isinstance(jd, (int, float)):
        raise TypeError("JD must be a number.")

    jd = float(jd)
    
    if jd < 2299160.5:
        raise ValueError("Dates before 1582-10-15 are not supported.")
    
    jd_shifted = jd + 0.5
    z = int(jd_shifted)
    f = jd_shifted - z

    if z < 2299161:
        a = z
    else:
        alpha = int((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - alpha // 4

    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)

    day_float = b - d - int(30.6001 * e) + f
    day = int(day_float)
    frac_day = day_float - day

    if e < 14:
        month = e - 1
    else:
        month = e - 13

    if month > 2:
        year = c - 4716
    else:
        year = c - 4715

    total_minutes = round(frac_day * 1440)

    if total_minutes >= 1440:
        total_minutes = 0
        base = datetime.datetime(year, month, day) + datetime.timedelta(days=1)
        year = base.year
        month = base.month
        day = base.day

    hour = total_minutes // 60
    minute = total_minutes % 60

    dt = datetime.datetime(year, month, day, hour, minute)

    return dt.strftime("%Y-%m-%dT%H:%M")