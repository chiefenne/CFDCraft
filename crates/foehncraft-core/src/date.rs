//! Calendar dates from Unix time, for the dates that reports and run
//! records carry (D84: taken from `SOURCE_DATE_EPOCH` when it is set, so
//! that two runs give the same files). Civil-from-days algorithm of
//! H. Hinnant, *chrono-compatible low-level date algorithms*.

/// Year, month, day, hour, minute, second (UTC) of `secs` seconds after
/// 1970-01-01T00:00:00Z.
pub fn civil(secs: i64) -> [i64; 6] {
    let z = secs.div_euclid(86_400) + 719_468;
    let s = secs.rem_euclid(86_400);
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    [y, m, d, s / 3600, s % 3600 / 60, s % 60]
}

#[cfg(test)]
mod tests {
    use super::civil;

    #[test]
    fn known_dates() {
        assert_eq!(civil(0), [1970, 1, 1, 0, 0, 0]);
        assert_eq!(civil(951_782_400), [2000, 2, 29, 0, 0, 0]);
        assert_eq!(civil(1_794_787_199), [2026, 11, 15, 23, 59, 59]);
        assert_eq!(civil(-1), [1969, 12, 31, 23, 59, 59]);
    }
}
