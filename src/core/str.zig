pub inline fn str_next_line(index: usize, slice: []const u8) usize {
    var idx = index;
    while (idx < slice.len) : (idx += 1) {
        if (slice[idx] == '\n') {
            idx += 1;
            return idx;
        } else if (slice[idx] == '\r') {
            idx += 1;
            if (idx < slice.len and slice[idx] == '\n') {
                idx += 1;
            }
            return idx;
        }
    }
    return slice.len;
}

pub inline fn str_r_trim(slice: []const u8) []const u8 {
    var result = slice;
    while (result.len > 0) {
        if (result[result.len - 1] == '\r' or result[result.len - 1] == '\n') {
            result = result[0..(result.len - 1)];
        }
        break;
    }
    return result;
}

pub inline fn is_carriage_return_char(c: u8) bool {
    return c == '\r' or c == '\n';
}

pub inline fn is_empty_char(c: u8) bool {
    return c == '\r' or c == '\n' or c == ' ';
}


