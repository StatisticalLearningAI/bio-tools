const std = @import("std");
const testing = std.testing;

const writer = @import("writer.zig");
const reader = @import("reader.zig");

pub const FqSlice = struct {
    name: [] const u8,
    qual: [] const u8,
    seq: [] const u8,
};

pub const FqArrayListRecord = struct {
    const Self = @This();
    const FqArrayList = std.ArrayListUnmanaged(u8); 
    
    allocator: std.mem.Allocator,
    name: FqArrayList,
    qual: FqArrayList,
    seq: FqArrayList,

    pub fn clear(self: *FqArrayListRecord) void {
        self.name.resize(0) catch {};
        self.qual.resize(0) catch {};
        self.seq.resize(0) catch {};
    }

    pub const FqWrite = FqWritable(
        *Self,
        resetRecord,
        FqArrayList.Writer,
        FqArrayList.Writer,
        FqArrayList.Writer
    );

    pub fn deinit(self: *FqArrayListRecord) void {
        self.name.deinit(self.allocator);
        self.qual.deinit(self.allocator);
        self.seq.deinit(self.allocator);
    }

    pub fn init(allocator: std.mem.Allocator) FqArrayListRecord {
        return .{
            .allocator = allocator,
            .name = .{},
            .qual = .{},
            .seq = .{},
        };
    }

    // writer to fill out the record
    pub fn writer(self: *Self) FqWrite  {
        return .{ 
            .context = self, 
            .name = self.name.writer(self.allocator),
            .qual = self.qual.writer(self.allocator),
            .seq = self.seq.writer(self.allocator),
        };
    }

    pub fn toSlice(self: *Self) FqSlice {
        return .{
            .name = self.name.items,
            .qual = self.qual.items,
            .seq = self.seq.items,
        };
    }

    pub fn resetRecord(self: *Self) void {
       self.name.shrinkRetainingCapacity(0); 
       self.qual.shrinkRetainingCapacity(0); 
       self.seq.shrinkRetainingCapacity(0); 
    }
};

pub fn FqWritable(comptime Context: type,
    comptime resetFn: fn(context: Context) void, // called to reset the record
    comptime NameWriter: type,
    comptime SeqWriter: type,
    comptime QualWriter: type,
) type {

    return struct {
        const Self = @This();
        
        context: Context,
        name: NameWriter,
        seq: SeqWriter,
        qual: QualWriter,
        
        pub fn resetRecord(self: Self) void {
           return resetFn(self.context); 
        }
    };
}

pub fn createReader(read: anytype) reader.FqReader(@TypeOf(read)) {
    return .{ .reader = read, .start_read_name = false};
}

pub fn createWriter(write: anytype) writer.FqWriter(@TypeOf(write)) {
    return .{.writer = write};
}


test {
    _ = reader;
    _ = writer;
}
