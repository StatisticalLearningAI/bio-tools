const std = @import("std");
const testing = std.testing;

const writer = @import("writer.zig");
const reader = @import("reader.zig");

pub const FqSlice = struct {
    name: []const u8,
    qual: []const u8,
    seq: []const u8,
};

// reads FQ data
//  Note: the data store in the result is temporary
//
pub fn fqReader(allocator: std.mem.Allocator ,stream: anytype, comptime option: reader.Option) reader.FqReader(@TypeOf(stream), option) {
    return .{
        .stream = stream,
        .stage = std.ArrayList(u8).init(allocator)
    };
}

pub fn fqWriter( stream: anytype, comptime option: writer.Option) writer.FqWriter(@TypeOf(stream.*), option) {
    return .{ 
        .stream = stream, 
    };
}

test {
    _ = reader;
    _ = writer;
}
