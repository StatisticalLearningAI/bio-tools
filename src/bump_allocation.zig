const std = @import("std");

pub const Alloc = struct {
    base_ptr: [*]u8,
    used: usize
};

// a temporary allocator that frees 
pub fn BumpAllocation(chunk_size: usize, alignment: usize) type {
    _ = chunk_size;
    _ = alignment;
    return struct {
        allocator: std.mem.Allocator,
        bump: ?Alloc = null,

        pub fn requestChunk(size: usize) ?[*]u8 {
            _ = size;
            return null;
        }

    };
}
