const std = @import("std");
const seq = @import("../seq.zig"); 
const testing = std.testing;

pub fn FqWriter(comptime Writer: type) type {
    return struct {
        const Self = @This();
        
        writer: Writer,

        pub fn write_record(self: *Self, slice: seq.fq.FqSlice) Writer.Error!void {
            try self.writer.writeByte('@');
            try self.writer.writeAll(slice.name);
            try self.writer.writeByte('\n');
            
            try self.writer.writeAll(slice.seq);
            try self.writer.writeByte('\n');
            
            try self.writer.writeByte('+');
            try self.writer.writeByte('\n');

            try self.writer.writeAll(slice.qual);
            try self.writer.writeByte('\n');
        }
    };
}



test "fq read test" {
    var buf = std.ArrayList(u8).init(std.testing.allocator);
    defer buf.deinit();
    var writer = seq.fq.createWriter(buf.writer());
    try writer.write_record(.{
        .name = "HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG",
        .seq = "TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGTATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT",
        .qual = "HIHIIIIIHIIHGHHIHHIIIIIIIIIIIIIIIHHIIIIIHHIHIIIIIGIHIIIIHHHHHHGHIHIIIIIIIIIIIGHIIIIIGHIIIIHIIHIHHIIIIHIHHIIIIIIIGIIIIIIIHIIIIIGHIIIIHIIIH?DGHEEGHIIIIIIIIIIIHIIHIIIHHIIHIHHIHCHHIIHGIHHHHHHH<GG?B@EHDE-BEHHHII5B@GHHF?CGEHHHDHIHIIH",
    });
    try testing.expectEqualStrings(buf.items,@embedFile("./test/t2.fq"));
}
