const std = @import("std");
const testing = std.testing;

pub const fai = @import("./fai.zig"); 
pub const fa = @import("./fa.zig");

pub const SeqContainer = @import("container.zig").SeqContianer;

pub const faiReader = fai.faiReaderIterator;

test {
    _ = fai;
    _ = fa;
}

