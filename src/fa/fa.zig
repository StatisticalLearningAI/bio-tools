const std = @import("std");
const testing = std.testing;

pub const FaSlice = struct {
    name: [] const u8,
    seq: [] const u8
};

pub const FaArrayList = struct {
    const Self = @This();

    allocator: std.mem.Allocator,
    name: std.ArrayListUnmanaged(u8) = . {},
    seq: std.ArrayListUnmanaged(u8)= .{},
};
