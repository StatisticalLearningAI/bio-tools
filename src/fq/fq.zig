const std = @import("std");
const testing = std.testing;

const writer = @import("writer.zig");
const reader = @import("reader.zig");

pub const FqSlice = struct {
    name: []const u8,
    qual: []const u8,
    seq: []const u8,
};

pub const FqArrayList = struct {
    const Self = @This();

    allocator: std.mem.Allocator,
    name: std.ArrayListUnmanaged(u8) = .{},
    qual: std.ArrayListUnmanaged(u8) = .{},
    seq: std.ArrayListUnmanaged(u8) = .{},

    pub const Error = error{OutOfMemory};


    pub fn deinit(self: *FqArrayList) void {
        self.name.deinit(self.allocator);
        self.qual.deinit(self.allocator);
        self.seq.deinit(self.allocator);
    }

    pub fn init(allocator: std.mem.Allocator) FqArrayList {
        return .{
            .allocator = allocator,
            .name = .{},
            .qual = .{},
            .seq = .{},
        };
    }

    pub fn toSlice(self: *Self) FqSlice {
        return .{
            .name = self.name.items,
            .qual = self.qual.items,
            .seq = self.seq.items,
        };
    }

    pub fn reset(self: *Self) void {
        self.name.shrinkRetainingCapacity(0);
        self.qual.shrinkRetainingCapacity(0);
        self.seq.shrinkRetainingCapacity(0);
    }

    pub fn nameWriter(self: *Self)  std.ArrayListUnmanaged(u8).Writer{
        return self.name.writer(self.allocator); 
    }
    
    pub fn seqWriter(self: *Self)  std.ArrayListUnmanaged(u8).Writer {
        return self.seq.writer(self.allocator); 
    }

    pub fn qualWriter(self: *Self) std.ArrayListUnmanaged(u8).Writer {
        return self.qual.writer(self.allocator); 
    }
};

pub fn createReader(read: anytype, comptime options: reader.FqReadOptions) reader.FqReadIterator(@TypeOf(read), options) {
    return .{ .reader = read };
}

pub fn createWriter(write: anytype, comptime options: writer.FqWriteOptions) writer.FqWriter(@TypeOf(write), options) {
    return .{ .writer = write };
}

test {
    _ = reader;
    _ = writer;
}
