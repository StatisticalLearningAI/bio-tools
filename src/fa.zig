const seq = @import("./seq.zig");
const std = @import("std");
const testing = std.testing;

pub const FqWriteOption = struct {
    lineEnding: enum { crlf, lf } = .lf,
};

pub fn FqWriter(comptime Stream: type) type {
    return struct {
        const Self = @This();

        stream: *Stream,
        options: FqWriteOption,
        pub const WriterError = Stream.Writer.Error;

        pub fn write_record(self: *Self, slice: struct {
            name: []const u8,
            qual: []const u8,
            seq: []const u8,
        }) Stream.Writer.Error!void {
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

pub const SectionType = enum { name, qual, seq };
// these are read in order with fastq
pub const Section = union(SectionType) {
    name: struct { format: seq.SeqContianer, data: []const u8 },
    qual: struct { line: u32, chunk: []const u8 },
    seq: struct { line: u32, chunk: []const u8 },
};

pub fn Writer(
    comptime Context: type,
    comptime writeFn: fn (context: Context, chunk: Section) anyerror!void,
) type {
    return struct {
        const Self = @This();
        context: Context,
        pub fn write(self: Self, chunk: Section) anyerror!void {
            try writeFn(self.context, chunk);
        }
    };
}

pub const FaRecord = struct {
    const Self = @This();
    buf: std.ArrayList(u8),
    format: seq.SeqContianer = .unknown,
    name_len: usize = 0,
    qual_len: usize = 0,
    seq_len: usize = 0,

    pub const FastaWriter = seq.fa.Writer(*Self, consumerHandler);

    pub fn deinit(self: *FaRecord) void {
        self.buf.deinit();
    }

    pub fn fastaWriter(self: *Self) FastaWriter { 
        return .{ .context = self };
    }

    fn consumerHandler(self: *Self, chunk: Section) anyerror!void {
        switch (chunk) {
            .name => |c| {
                try self.buf.appendSlice(c.data);
                self.name_len += c.data.len;
                self.format = c.format;
            },
            .qual => |c| {
                try self.buf.appendSlice(c.chunk);
                self.qual_len += c.chunk.len;
            },
            .seq => |c| {
                try self.buf.appendSlice(c.chunk);
                self.seq_len += c.chunk.len;
            },
        }
    }

    fn name(self: *Self) []const u8 {
        return self.buf.items[0..self.name_len]; 
    }
    
    fn sequence(self: *Self) []const u8 {
        return self.buf.items[self.name_len .. self.name_len + self.seq_len];
    }

    fn quality(self: *Self) []const u8 {
        return self.buf.items[self.name_len + self.seq_len .. self.qual_len + self.seq_len + self.name_len]; 
    }

    fn reset(self: *Self) void {
        self.buf.shrinkRetainingCapacity(0);
        self.name_len = 0;
        self.seq_len = 0;
        self.qual_len = 0;
    }
};

pub fn FaReadIterator(comptime Stream: type) type {
    return struct {
        const Self = @This();

        stream: Stream,
        
        state: enum { start, name, seq, qual } = .start,
        format: seq.SeqContianer = .fastq,

        pub const Option = struct { 
            reserve_len: u32 = 2048
            
        };
        pub fn next(self: *Self, consumer: anytype, comptime option: Option) (anyerror || error{BufferOverflow})!bool {
            var stage_buffer: [option.reserve_len]u8 = undefined;
            var stage = std.io.fixedBufferStream(stage_buffer[0..]);

            var reader = self.stream.reader();

            var seq_line_count: u32 = 0;
            var qual_line_count: u32 = 0;

            var seq_number_bases: usize = 0;
            var qual_number_bases: usize = 0;

            var bytes_per_line: usize = 0;
            var bytes_per_line_expected: ?usize  = null;
            
            var bases_per_line: usize = 0;
            var bases_per_line_expected: ?usize = null;

            while (true) {
                switch (self.state) {
                    .start => {
                        const byte: u8 = reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => return false, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        self.state = switch (byte) {
                            '@' => l: {
                                self.format = .fastq;
                                break :l .name;
                            },
                            '>' => l: {
                                self.format = .fasta;
                                break :l .name;
                            },
                            else => return error.RecordParseError,
                        };
                    },
                    .name => {
                        reader.streamUntilDelimiter(stage.writer(), '\n', option.reserve_len) catch |err| switch (err) {
                            error.EndOfStream => return false,
                            else => |e| return e,
                        };
                        try consumer.write(.{ .name = .{ .format = self.format, .data = std.mem.trim(u8, stage.getWritten(), &std.ascii.whitespace) } });
                        stage.reset();
                        self.state = .seq;
                    },
                    .seq => blk: {
                        while (true) {
                            switch (try reader.readByte()) {
                                '+' => {
                                    self.state = .qual;
                                    try reader.skipUntilDelimiterOrEof('\n');
                                    bases_per_line = 0; 
                                    bytes_per_line = 0;
                                    break :blk;
                                },
                                '>' => {
                                    if (self.format == .fastq)
                                        return error.ParseError;
                                    self.state = .name; // return back to the name state
                                    return true;
                                },
                                else => |b| {
                                    if(seq_line_count > 0) {
                                        // we're onto the next line need to check line_counts
                                        if(bytes_per_line_expected) |count| {
                                            if(bytes_per_line != count) return error.RecordParseError;
                                        } else bytes_per_line_expected = bytes_per_line;
                                        if(bases_per_line_expected) |count| {
                                            if(bases_per_line != count) return error.RecordParseError;
                                        } else bases_per_line_expected = bases_per_line;
                                    }
                                        
                                    bases_per_line = 0; 
                                    bytes_per_line = 0;

                                    const app: [1]u8 = .{b};
                                    _ = try stage.write(&app);
                                },
                            }
                            same_line_chunk: while(true)  {
                                seq_line_count += 1;
                                reader.streamUntilDelimiter(stage.writer(), '\n', option.reserve_len) catch |err| switch (err) {
                                    error.EndOfStream => return false,
                                    error.StreamTooLong => {
                                        seq_line_count -= 1;
                                        continue :same_line_chunk;
                                    },
                                    error.NoSpaceLeft => unreachable,
                                    else => |e| return e,
                                };
                                const byte_buf = stage.getWritten();
                                const bases_buf = std.mem.trim(u8, byte_buf, &std.ascii.whitespace);
                              
                                bytes_per_line += byte_buf.len;
                                bases_per_line += bases_buf.len;

                                try consumer.write(.{ .seq = .{ .line = seq_line_count, .chunk = bases_buf } });
                                
                                seq_number_bases += bases_buf.len;
                                stage.reset();
                                break :same_line_chunk;
                            }
                        }
                    },
                    .qual => {

                        if(qual_line_count > 0) {
                            if(bytes_per_line_expected) |count| {
                                if(bytes_per_line != count) return error.RecordParseError;
                            } else return error.RecordParseError;
                            if(bases_per_line_expected) |count| {
                                if(bases_per_line != count) return error.RecordParseError;
                            } else return error.RecordParseError; 
                            
                            bases_per_line = 0; 
                            bytes_per_line = 0;
                        }

                        same_line_chunk: while(true)  {
                            qual_line_count += 1;
                            reader.streamUntilDelimiter(stage.writer(), '\n', option.reserve_len) catch |err| switch (err) {
                                error.EndOfStream => return false,
                                error.StreamTooLong => {
                                    qual_line_count -= 1;
                                    continue :same_line_chunk;
                                },
                                error.NoSpaceLeft => unreachable,
                                else => |e| return e,
                            };
                            const byte_buf = stage.getWritten();
                            const qual_buf = std.mem.trim(u8, byte_buf, &std.ascii.whitespace);
                            
                            bytes_per_line += byte_buf.len;
                            bases_per_line += qual_buf.len;

                            try consumer.write(.{ .qual = .{ .line = qual_line_count, .chunk = qual_buf } });
                            qual_number_bases += qual_buf.len;
                            stage.reset();
                            break :same_line_chunk;
                        }
                       
                        if (qual_number_bases == seq_number_bases) { 
                            if(seq_line_count != qual_line_count) return error.RecordParseError; 
                            self.state = .start;
                            return true;
                        } else if (qual_number_bases > seq_number_bases) {
                            return error.RecordParseError;
                        }
                    },
                }
            }
            return false;
        }
    };
}

pub fn fastaReadIterator( stream: anytype) FaReadIterator(@TypeOf(stream)) {
    return .{ .stream = stream };
}
//pub fn readUsingIndex(self: *Self, index: anytype) void {
//
//
//}

pub fn FqReader(comptime Stream: type) type {
    return struct {
        stream: Stream,
        stage: std.ArrayList(u8),
        const Self = @This();

        pub const ParseErrors = error{RecordParseError};
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

        pub const ParseOptions = struct {};

        pub fn next(self: *Self, comptime options: ParseOptions) (Stream.Reader.Error || ParseErrors || error{OutOfMemory})!?struct {
            length: u32, // Total length of this reference sequence, in bases
            offset: usize, // Offset in the FASTA/FASTQ file of this sequence's first base
            line_base: u32, // The number of bases on each line
            line_width: u32, // The number of bytes in each line, including the newline
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
                            '\n', '\r' => return error.RecordParseError,
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
                    .seq => {
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
                                if (expected_bases_per_line == null) expected_bases_per_line = bases_per_line;
                                if (expected_bytes_per_line == null) expected_bytes_per_line = bytes_per_line;
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
                                if (state == .qual_start)
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
                return error.RecordParseError;

            return .{ .length = seq_len, .offset = seq_read_offset, .line_base = bases_per_line, .line_width = bytes_per_line, .qual_offset = qual_read_offset, 
                .name = self.stage.items[0..name_len], 
                .seq = self.stage.items[name_len .. name_len + seq_len], 
                .qual = self.stage.items[name_len + seq_len .. qual_len + seq_len + name_len] };
        }
    };
}

pub fn fqReader(allocator: std.mem.Allocator, stream: anytype) FqReader(@TypeOf(stream)) {
    return .{ .stream = stream, .stage = std.ArrayList(u8).init(allocator) };
}

//test "parse fai index and seek fq" {
//    //const indexFile = @embedFile("./test/t3.fq.fai");
//    const fq = @embedFile("./test/t3.fq");
//    var fqStream = std.io.fixedBufferStream(fq);
//    var fqReadIter = seq.fq.fqReader(std.testing.allocator, fqStream);
//    defer fqReadIter.deinit();
//    const TestCase = struct { index: seq.fai.FaiEntry(.fastq, .{}), seq: []const u8, qual: []const u8 };
//
//    const test_cases =
//        [_]TestCase{
//        .{ .index = .{ .name = "FAKE0005_1", .length = 63, .offset = 85, .line_base = 63, .line_width = 64, .qual_offset = 151 }, .seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG", .qual = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~" },
//        .{ .index = .{ .name = "FAKE0006_4", .length = 63, .offset = 1590, .line_base = 63, .line_width = 64, .qual_offset = 1656 }, .seq = "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA", .qual = "~}|{zyxwvutsrqponmlkjihgfedcba`_^]\\[ZYXWVUTSRQPONMLKJIHGFEDCBA@" },
//    };
//
//    for (test_cases) |case| {
//        const record = try fqReadIter.readUsingIndex(case.index, .{});
//        try testing.expectEqualStrings(case.seq, record.seq);
//        try testing.expectEqualStrings(case.qual, record.qual);
//    }
//}

test "read fq records t1.fq" {
    const file = @embedFile("./test/t1.fq");
    var stream = std.io.fixedBufferStream(file);
    var iter = seq.fa.fastaReadIterator(stream);
    const TestCase = struct { 
        name: []const u8,
        sequence: []const u8, 
        quality: []const u8 
    };

    const test_cases =
        [_]TestCase{
        .{
            .name = "HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG",
            .sequence = "TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGTATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT",
            .quality = "HIHIIIIIHIIHGHHIHHIIIIIIIIIIIIIIIHHIIIIIHHIHIIIIIGIHIIIIHHHHHHGHIHIIIIIIIIIIIGHIIIIIGHIIIIHIIHIHHIIIIHIHHIIIIIIIGIIIIIIIHIIIIIGHIIIIHIIIH?DGHEEGHIIIIIIIIIIIHIIHIIIHHIIHIHHIHCHHIIHGIHHHHHHH<GG?B@EHDE-BEHHHII5B@GHHF?CGEHHHDHIHIIH"
        }, 
        .{

            .name = "HWI-D00523:240:HF3WGBCXX:1:1101:5586:3020 1:N:0:CTGTAG",
            .sequence = "TGGGGAATATTGGGCAATGGGCGGAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCCCTCGGGTTGTAAACTTCTTTTCTATAGGACGAAGAAGTGACGGTACTATAGGAATAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCGAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGAGAGCAAGTCAGATGTGA",
            .quality = "EHEHGIIIHIGGHGHFEHHEHGCHHGGHIIIGHHIFHHGHHIEHIIIIGHHHHIIGIHGHGGHHHHHCHHHICCHHHHH@HHIIIGCEHGHHGHCHHGDGCGCCEHEEHGIIGHHGHHHIGGFCFHHIHIGIIIHGGHFHIIFEFHIIHIGDHCHFHHGCHCE?GHIIH<C?GHHHHIGFDEHHHHEC88@<@@EHHHIHDH-@HHCDHHDDEHH6@F6@6@@EH@@"
        }, 
        .{

            .name = "HWI-D00523:240:HF3WGBCXX:1:1101:2860:2149 1:N:0:CTGTAG",
            .sequence = "TAGGGAATATTGCTCAATGGGGGAAACCCTGAAGCAGCAACGCCGCGTGGAGGATGAAGGTTTTAGGATTGTAAACTCCTTTTGTGAGAGAAGATTATGACGGTATCTCACGAATAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCGAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGTAGGCGGGACTGCAAGTTGGGTGTCAAA",
            .quality = "HGHHGHIIHIIIIIIHHHIIHIIIIIHGHCHIHIIHIIHIIIIIIIIIHIGHIEHIHIIG<FEHHHIHHIIIIIIHIFFHHHHIIIIIHHHHGHFEHHIHHHIHEHHFHHHIHIIIIIIIHIHDHHHHHEHHIIGIIHHIIGHHHIIGDDAGGHHFHHHIHICHHHGH,GHHHHGCEHEG?@6@G?-@>HHHHHHHDEH<@H-@CDD>:E?@GHEF-@E:@H+@-@@"
        }
    };

    var record: FaRecord = .{ .buf = std.ArrayList(u8).init(std.testing.allocator) };
    defer record.deinit();
    for (test_cases) |case| {
        if(try iter.next(record.fastaWriter(), .{})) {
            try testing.expectEqualStrings(case.name, record.name());
            try testing.expectEqualStrings(case.sequence, record.sequence());
            try testing.expectEqualStrings(case.quality, record.quality());
        } else try testing.expect(false);
        record.reset();
    }

    // defer iter.deinit();
   // if (try iter.next(record.fastaWriter(), .{}) {
   //     try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG", result.name);
   //     try testing.expectEqualStrings("TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGTATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT", result.seq);
   // } else try testing.expect(false);

   // if (try fqReadIter.next(.{})) |result| {
   //     try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:5586:3020 1:N:0:CTGTAG", result.name);
   // } else try testing.expect(false);

   // if (try fqReadIter.next(.{})) |result| {
   //     try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:2860:2149 1:N:0:CTGTAG", result.name);
   // } else try testing.expect(false);
}

//test "writing fq to buffer" {
//    var buf = std.ArrayList(u8).init(std.testing.allocator);
//    defer buf.deinit();
//    var writer = seq.fq.fqWriter(&buf, .{});
//    try writer.write_record(.{
//        .name = "HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG",
//        .seq = "TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGTATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT",
//        .qual = "HIHIIIIIHIIHGHHIHHIIIIIIIIIIIIIIIHHIIIIIHHIHIIIIIGIHIIIIHHHHHHGHIHIIIIIIIIIIIGHIIIIIGHIIIIHIIHIHHIIIIHIHHIIIIIIIGIIIIIIIHIIIIIGHIIIIHIIIH?DGHEEGHIIIIIIIIIIIHIIHIIIHHIIHIHHIHCHHIIHGIHHHHHHH<GG?B@EHDE-BEHHHII5B@GHHF?CGEHHHDHIHIIH",
//    });
//    try testing.expectEqualStrings(buf.items, @embedFile("./test/t2.fq"));
//}
