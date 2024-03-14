const seq = @import("./seq.zig");
const std = @import("std");
const testing = std.testing;

pub fn FaReader(comptime Stream: type) type {
    return struct {
        stream: Stream,
        stage: std.ArrayList(u8),
        const Self = @This();
    };

}

pub fn faReader(allocator: std.mem.Allocator, stream: anytype) FaReader(@TypeOf(stream)) {
    return .{ .stream = stream, .stage = std.ArrayList(u8).init(allocator) };
}
