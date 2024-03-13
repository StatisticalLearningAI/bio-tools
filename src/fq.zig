const seq = @import("./seq.zig");
const std = @import("std");
const testing = std.testing;

pub const FqEntry = struct {
    name: []const u8,
    qual: []const u8,
    seq: []const u8,
};

pub const FqWriteOption = struct { 
    repeatNameForQual: bool = false, 
};

pub fn FqWriter(comptime Stream: type, comptime options: FqWriteOption) type {

    return struct {
        const Self = @This();
        stream: *Stream,
        pub const WriterError = Stream.Writer.Error;
        pub fn write_record(self: *Self, slice: seq.fq.FqEntry) Stream.Writer.Error!void {
            var fq = if(@hasField(Stream, "fq")) self.stream.fq.writer() else self.stream.writer();

            try fq.writeByte('@');
            try fq.writeAll(slice.name);
            try fq.writeByte('\n');

            try fq.writeAll(slice.seq);
            try fq.writeByte('\n');

            try fq.writeByte('+');
            if (options.repeatNameForQual) {
                try fq.writeAll(slice.name);
            }
            try fq.writeByte('\n');

            try fq.writeAll(slice.qual);
            try fq.writeByte('\n');
        }
    };
}

pub fn fqWriter( stream: anytype, comptime option: FqWriteOption) FqWriter(@TypeOf(stream.*), option) {
    return .{ 
        .stream = stream, 
    };
}


pub const FqReaderOption = struct {
    skipQual: bool = false, // skip qualifiers 
    //lineCountOfQualAndSeqMatch: bool = false, // ensure the number of lines matches for qual and seq
};

pub fn FqReader(comptime Stream: type, comptime options: FqReaderOption) type {
    return struct {
        stream: Stream,
        stage: std.ArrayList(u8),  
        const Self = @This();

        pub const Container = seq.SeqContianer.fastq;

        pub const FqIterator = struct {
            seq_offset: usize, // offset in the reader
            sequence_len: usize, // the number of bases in the sequence
            qual_offset: u64,
            entry: seq.fq.FqEntry
        }; 

        pub fn deinit(self: *Self) void {
            self.stage.deinit();
        }

        // seek to the record using FaiEntry
        pub fn readUsingFai(self: *Self, fai: anytype) (Stream.SeekableStream.Error)!void{
            _ = fai;
            var seek = self.stream.seekableStream();
            _ = seek;
        }

        pub fn next(self: *Self) (Stream.Reader.Error || error{ OutOfMemory, MismatchedQualAndSeq, MalformedIdentifier })!?FqIterator {
            self.stage.shrinkRetainingCapacity(0);
            var seq_read_offset: usize = 0;
            var qual_read_offset: usize = 0;

            var reader = self.stream.reader();

            var count_seq: u32 = 0;
            var count_qul: u32 = 0;
            var state: enum { start, name, seq, seq_start, qual, qual_start } = .start;
            
            var name_len: usize = 0;
            var seq_len: usize  = 0;
            var qual_len: usize  = 0;

            finished: while (true) {
                switch (state) {
                    .start => {
                        const byte: u8 = reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => return null, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        state = switch (byte) {
                            '@' => .name,
                            else => return error.MalformedIdentifier,
                        };
                    },
                    .name => {
                        const byte: u8 = reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => return null, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        state = switch (byte) {
                            '\r' => state, // ignore \r characters
                            '\n' => .seq_start, // end of the line
                            else => l: {
                                try self.stage.append(byte);
                                name_len += 1;
                                break :l state;
                            },
                        };
                    },
                    .seq_start => {
                        const byte: u8 = reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => continue :finished, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        state = switch (byte) {
                            '+' => .qual_start,
                            else => l: {
                                seq_read_offset = try self.stream.seekableStream().getPos();
                                count_seq += 1;
                                try self.stage.append(byte);
                                seq_len += 1;
                                break :l .seq;
                            },
                        };
                    },
                    .seq => {
                        const byte: u8 = reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => return null, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        state = switch (byte) {
                            '\r' => state,
                            '\n' => .seq_start,
                            else => l: {
                                count_seq += 1;
                                try self.stage.append(byte);
                                seq_len += 1;
                                break :l state;
                            },
                        };
                    },
                    .qual, .qual_start => {
                        const byte: u8 = reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => return null, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        state = switch (byte) {
                            '\r' => .qual,
                            '\n' => l: {
                                // '@' and '+' can show up in qualifiers
                                //   the number of sequence characters has to equal qualifiers
                                if (count_qul == count_seq)
                                    break :finished;
                                if (count_qul > count_seq)
                                    return error.MismatchedQualAndSeq;
                                break :l .qual;
                            },
                            else => l: {
                                count_qul += 1;
                                if(state == .qual_start) {
                                    qual_read_offset = try self.stream.seekableStream().getPos(); 
                                }
                                //try qual_writer.writeByte(byte);
                                if (!comptime options.skipQual) {
                                    try self.stage.append(byte);
                                    qual_len += 1; 
                                }
                                break :l .qual;
                            },
                        };
                    },
                }
            }
            return .{ 
                .seq_offset = seq_read_offset, 
                .qual_offset = qual_read_offset,
                .sequence_len = count_seq, 
                .entry = .{
                    .name = self.stage.items[0..name_len],
                    .seq = self.stage.items[name_len..name_len + seq_len],
                    .qual = self.stage.items[name_len + seq_len..qual_len + seq_len + name_len]
                }
            };
        }
    };
}

pub fn fqReader(allocator: std.mem.Allocator ,stream: anytype, comptime option: FqReaderOption) FqReader(@TypeOf(stream), option) {
    return .{
        .stream = stream,
        .stage = std.ArrayList(u8).init(allocator)
    };
}

test "test read fq records" {
    const file = @embedFile("./test/t1.fq");
    var stream = std.io.fixedBufferStream(file);
    var fqReadIter = seq.fq.fqReader(std.testing.allocator,stream, .{});
    defer fqReadIter.deinit();
    if (try fqReadIter.next()) |result| {
        try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG", result.entry.name);
        try testing.expectEqualStrings("TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGTATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT", result.entry.seq);
    } else try testing.expect(false);

    if (try fqReadIter.next()) |result| {
        try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:5586:3020 1:N:0:CTGTAG", result.entry.name);
    } else try testing.expect(false);

    if (try fqReadIter.next()) |result| {
        try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:2860:2149 1:N:0:CTGTAG", result.entry.name);
    } else try testing.expect(false);
}

test "fq writer test" {
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
