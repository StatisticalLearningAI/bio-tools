const seq = @import("./seq.zig");
const std = @import("std");
const testing = std.testing;

pub const FqWriteOption = struct {
    lineEnding: enum {
       crlf,
       lf 
    } = .lf,
};

pub fn FqWriter(comptime Stream: type) type {
    return struct {
        const Self = @This();

        stream: *Stream,
        options: FqWriteOption,
        pub const WriterError = Stream.Writer.Error;

        pub fn write_record(self: *Self, slice: anytype) Stream.Writer.Error!void {
            var fq = self.stream.writer();

            try fq.writeByte('@');
            try fq.writeAll(slice.name);
            try fq.writeByte('\n');

            try fq.writeAll(slice.seq);
            try fq.writeByte('\n');

            try fq.writeByte('+');
            try fq.writeByte('\n');

            try fq.writeAll(slice.qual);
            try fq.writeByte('\n');
        }
    };
}

pub fn fqWriter(stream: anytype, option: FqWriteOption) FqWriter(@TypeOf(stream.*)) {
    return .{
        .options = option,
        .stream = stream,
    };
}

pub fn FqReader(comptime Stream: type) type {
    return struct {
        stream: Stream,
        stage: std.ArrayList(u8),
        const Self = @This();

        pub const ParseErrors = error{ RecordParseError };
        pub const Container = seq.SeqContianer.fastq;

        pub fn deinit(self: *Self) void {
            self.stage.deinit();
        }

        // optional parameters for readining index
        pub const ReadUsingIndexOptions = struct { skip_qual: bool = false };
        // seek to the record using FaiEntry
        pub fn readUsingIndex(self: *Self, fai: anytype, comptime options: ReadUsingIndexOptions) (Stream.SeekError || error{ OutOfMemory, EndOfStream })!struct {
            name: []const u8,
            qual: []const u8,
            seq: []const u8,
        } {
            self.stage.shrinkRetainingCapacity(0);
            var seek = self.stream.seekableStream();
            var reader = self.stream.reader();

            {
                try seek.seekTo(fai.offset); // seek to the offset of the fai
                var remaining_bases: usize = fai.length;
                while (remaining_bases > 0) {
                    const bases_to_read = @min(fai.line_base, remaining_bases);
                    try self.stage.ensureTotalCapacity(self.stage.items.len + bases_to_read);
                    const bytes_read = try reader.read(self.stage.unusedCapacitySlice()[0..bases_to_read]);
                    if (bytes_read < bases_to_read)
                        return error.EndOfStream;
                    self.stage.items.len += bytes_read;
                    remaining_bases -= bytes_read;
                    try reader.skipBytes(fai.line_width - fai.line_base, .{ .buf_size = 16 });
                }
            }
            if (!options.skip_qual) {
                try seek.seekTo(fai.qual_offset); // seek to the offset of the fai
                var remaining_qual: usize = fai.length;
                while (remaining_qual > 0) {
                    const bases_to_read = @min(fai.line_base, remaining_qual);
                    try self.stage.ensureTotalCapacity(self.stage.items.len + bases_to_read);
                    const bytes_read = try reader.read(self.stage.unusedCapacitySlice()[0..bases_to_read]);
                    if (bytes_read < bases_to_read)
                        return error.EndOfStream;
                    self.stage.items.len += bytes_read;
                    remaining_qual -= bytes_read;
                    try reader.skipBytes(fai.line_width - fai.line_base, .{ .buf_size = 16 });
                }
            }
            return .{
                .name = fai.name,
                .seq = self.stage.items[0..fai.length],
                .qual = if (options.skip_qual) .{} else self.stage.items[fai.length .. fai.length + fai.length],
            };
        }

        pub const ParseOptions = struct {

        };

        pub fn next(self: *Self, comptime options: ParseOptions) (Stream.Reader.Error || ParseErrors || error{OutOfMemory})!?struct {
            length: u32,      // Total length of this reference sequence, in bases
            offset: usize,      // Offset in the FASTA/FASTQ file of this sequence's first base
            line_base: u32,   // The number of bases on each line
            line_width: u32,  // The number of bytes in each line, including the newline
            qual_offset: usize, // Offset of sequence's first quality within the FASTQ file

            name: []const u8,
            qual: []const u8,
            seq: []const u8,
        } {
            _ = options;
            self.stage.shrinkRetainingCapacity(0);
            var seq_read_offset: usize = 0;
            var qual_read_offset: usize = 0;

            var reader = self.stream.reader();

            var state: enum { start, name, seq, seq_begin, qual, qual_start } = .start;

            var name_len: usize = 0;
            var seq_len: u32 = 0;
            var qual_len: u64 = 0;

            var bases_per_line: u32 = 0;
            var expected_bases_per_line: ?u32 = null;

            var bytes_per_line: u32 = 0;
            var expected_bytes_per_line: ?u32 = null;

            var seq_num_lines: usize = 0;
            var qual_num_lines: usize = 0;

            finished: while (true) {
                switch (state) {
                    .start => {
                        const byte: u8 = reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => return null, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        state = switch (byte) {
                            '@' => .name,
                            else => return error.RecordParseError,
                        };
                    },
                    .name => {
                        const byte: u8 = reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => return null, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        state = switch (byte) {
                            '\r' => state, // ignore \r characters
                            '\n' => .seq_begin, // end of the line
                            else => l: {
                                try self.stage.append(byte);
                                name_len += 1;
                                break :l state;
                            },
                        };
                    },
                    .seq_begin => {
                        const byte: u8 = reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => '\n', // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        state = switch (byte) {
                            '+' => .qual_start,
                            '\n', '\r' =>  return error.RecordParseError,
                            else => l: {
                                if (expected_bytes_per_line) |e| 
                                    if (e != bytes_per_line) return error.RecordParseError;
                                 
                                if (expected_bases_per_line) |e| 
                                    if (e != bases_per_line) return error.RecordParseError;

                                bases_per_line = 1;
                                bytes_per_line = 1;
                                seq_len += 1;
                                try self.stage.append(byte);
                                
                                break :l .seq;
                            },
                        };
                    },
                    .seq =>  {

                        const byte: u8 = reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => return null, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        bytes_per_line += 1;
                        state = switch (byte) {
                            '+' => .qual_start,
                            '\r' => state,
                            '\n' => l: {
                                seq_num_lines += 1;
                                if (expected_bases_per_line == null ) expected_bases_per_line = bases_per_line;
                                if (expected_bytes_per_line == null ) expected_bytes_per_line = bytes_per_line;
                                break :l .seq_begin;
                            },
                            else => l: {
                                bases_per_line += 1;
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
                                if (state == .qual) qual_num_lines += 1;
                                // '@' and '+' can show up in qualifiers
                                //   the number of sequence characters has to equal qualifiers
                                if (qual_len == seq_len) break :finished;
                                if (qual_len > seq_len) return error.RecordParseError;
                                break :l .qual;
                            },
                            else => l: {
                                if (state == .qual_start )
                                    qual_read_offset = try self.stream.seekableStream().getPos();

                                try self.stage.append(byte);
                                qual_len += 1;
                                break :l .qual;
                            },
                        };
                    },
                }
            }

            if (qual_num_lines != seq_num_lines)
                return error.RecordParseError ;

            return .{ 
                .length = seq_len,
                .offset = seq_read_offset,
                .line_base = bases_per_line,
                .line_width = bytes_per_line,
                .qual_offset = qual_read_offset,

                .name = self.stage.items[0..name_len], 
                .seq = self.stage.items[name_len .. name_len + seq_len], 
                .qual = self.stage.items[name_len + seq_len .. qual_len + seq_len + name_len] 
            };
        }
    };
}

pub fn fqReader(allocator: std.mem.Allocator, stream: anytype) FqReader(@TypeOf(stream)) {
    return .{ .stream = stream, .stage = std.ArrayList(u8).init(allocator) };
}

test "parse fai index and seek fq" {
    //const indexFile = @embedFile("./test/t3.fq.fai");
    const fq = @embedFile("./test/t3.fq");
    var fqStream = std.io.fixedBufferStream(fq);
    var fqReadIter = seq.fq.fqReader(std.testing.allocator, fqStream);
    defer fqReadIter.deinit();
    const TestCase = struct { index: seq.fai.FaiEntry(.fastq, .{}), seq: []const u8, qual: []const u8 };

    const test_cases =
        [_]TestCase{
        .{ .index = .{ .name = "FAKE0005_1", .length = 63, .offset = 85, .line_base = 63, .line_width = 64, .qual_offset = 151 }, .seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG", .qual = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~" },
        .{ .index = .{ .name = "FAKE0006_4", .length = 63, .offset = 1590, .line_base = 63, .line_width = 64, .qual_offset = 1656 }, .seq = "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA", .qual = "~}|{zyxwvutsrqponmlkjihgfedcba`_^]\\[ZYXWVUTSRQPONMLKJIHGFEDCBA@" },
    };

    for (test_cases) |case| {
        const record = try fqReadIter.readUsingIndex(case.index, .{});
        try testing.expectEqualStrings(case.seq, record.seq);
        try testing.expectEqualStrings(case.qual, record.qual);
    }
}

test "read fq records t1.fq" {
    const file = @embedFile("./test/t1.fq");
    var stream = std.io.fixedBufferStream(file);
    var fqReadIter = seq.fq.fqReader(std.testing.allocator, stream);
    defer fqReadIter.deinit();
    if (try fqReadIter.next(.{})) |result| {
        try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG", result.name);
        try testing.expectEqualStrings("TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGTATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT", result.seq);
    } else try testing.expect(false);

    if (try fqReadIter.next(.{})) |result| {
        try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:5586:3020 1:N:0:CTGTAG", result.name);
    } else try testing.expect(false);

    if (try fqReadIter.next(.{})) |result| {
        try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:2860:2149 1:N:0:CTGTAG", result.name);
    } else try testing.expect(false);
}

test "writing fq to buffer" {
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
