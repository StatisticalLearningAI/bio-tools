const std = @import("std");
const testing = std.testing;

pub const str = @import("./core/str.zig");
pub const fq = @import("fq/fq.zig"); 
pub const fai = @import("fai/fai.zig"); 


test {
    _ = fq;
}
