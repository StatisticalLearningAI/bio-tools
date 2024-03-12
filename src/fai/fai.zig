const std = @import("std");

//Using an fai index file in conjunction with a FASTA/FASTQ 
pub const FaiSlice = struct {
    name: [] const u8,
    length: usize,
    offset: usize,
    line_base: usize,
    line_width: usize,
    qual_offset: usize
};

pub const FaiArrayList = struct {
    const Self = @This();
    allocator: std.mem.Allocator,
    
    name: std.ArrayListUnmanaged(u8),
    length: usize = 0,
    offset: usize = 0,
    line_base: usize = 0,
    line_width: usize = 0,
    qual_offset: usize = 0,

    pub fn init(allocator: std.mem.Allocator) FaiArrayList {
        return .{
            .allocator = allocator,
            .name = std.ArrayListUnmanaged(u8).init(allocator),
        };
    }

    pub fn deinit(self: *Self) void {
        self.name.deinit(self.allocator);
    }

    pub fn toSlice(self: *Self) FaiSlice {
        return .{
            .name = self.name.items,
            .length = self.length,
            .offset = self.offset,
            .line_base =  self.line_base,
            .line_width = self.line_width,
            .qual_offset = self.qual_offset 
        };
    }

    pub fn nameWriter(self: *Self) std.ArrayListUnmanaged(u8).Writer {
        return self.name.writer(self.allocator);
    }
    
    pub fn reset(self: *Self) void {
        self.name.shrinkRetainingCapacity(0);
    }

    pub fn set(self: *Self, len: usize, offset: usize, line_base: usize, line_width: usize, qual_offset: usize) void{
        self.length = len;
        self.offset = offset;
        self.line_base = line_base;
        self.line_width = line_width;
        self.qual_offset = qual_offset;
    } 
};

