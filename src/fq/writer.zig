const std = @import("std");
const seq = @import("../seq.zig");
const testing = std.testing;

//http://www.htslib.org/doc/faidx.html
pub fn FqStreamWithFai(comptime FQStream: type, comptime FQIStream: type) type {
    return struct {
        pub const Error = FQStream.Writer.Error || FQIStream.Writer.Error;
        fq: FQStream,
        fqI: FQIStream,
    };
}

pub const Option = struct { 
    addNameToQualOp: bool = false, 
};

pub fn FqWriter(comptime Stream: type, comptime options: Option) type {
    const isStreamWithFai = @hasField(Stream, "fq") and @hasField(Stream, "fqI");
    _ = isStreamWithFai;  

    return struct {
        const Self = @This();
        stream: *Stream,
        offset: usize = 0,

        pub const WriterError = (if(@hasField(Stream, "fq")) Stream.fq.Writer.Error else Stream.Writer.Error);

        pub fn write_record(self: *Self, slice: seq.fq.FqSlice) Stream.Writer.Error!void {
            var fq = if(@hasField(Stream, "fq")) self.stream.fq.writer() else self.stream.writer();

            try fq.writeByte('@');
            try fq.writeAll(slice.name); self.offset += slice.name.len;
            try fq.writeByte('\n');

            try fq.writeAll(slice.seq);
            try fq.writeByte('\n');

            try fq.writeByte('+');
            if (options.addNameToQualOp) {
                try fq.writeAll(slice.name);
            }
            try fq.writeByte('\n');

            try fq.writeAll(slice.qual);
            try fq.writeByte('\n');
        }
    };
}

test "fq read test" {
    var buf = std.ArrayList(u8).init(std.testing.allocator);
    defer buf.deinit();
    var writer = seq.fq.fqWriter(&buf, .{});
    try writer.write_record(.{
        .name = "HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG",
        .seq = "TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGTATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT",
        .qual = "HIHIIIIIHIIHGHHIHHIIIIIIIIIIIIIIIHHIIIIIHHIHIIIIIGIHIIIIHHHHHHGHIHIIIIIIIIIIIGHIIIIIGHIIIIHIIHIHHIIIIHIHHIIIIIIIGIIIIIIIHIIIIIGHIIIIHIIIH?DGHEEGHIIIIIIIIIIIHIIHIIIHHIIHIHHIHCHHIIHGIHHHHHHH<GG?B@EHDE-BEHHHII5B@GHHF?CGEHHHDHIHIIH",
    });
    try testing.expectEqualStrings(buf.items, @embedFile("./test/t2.fq"));
}
