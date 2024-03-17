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

pub const FaReadIterOptions = struct {
    reserve_len: usize = 4096,
};

pub fn FaReadIterator(comptime Stream: type, comptime options: FaReadIterOptions) type {
    return struct {
        const Self = @This();


        stream: Stream,
        buffer: [options.reserve_len]u8 = undefined,
        state: enum { start, name, seq, qual } = .start,
        format: seq.SeqContainer = .fastq,

        pub fn next(self: *Self, consumer: anytype) (@TypeOf(consumer).Error || Stream.ReadError || error{ BufferOverflow, RecordParseError })!bool {
            var stage = std.io.fixedBufferStream(&self.buffer);

            var reader = self.stream.reader();

            var seq_line_count: u32 = 0;
            var qual_line_count: u32 = 0;

            var seq_number_bases: usize = 0;
            var qual_number_bases: usize = 0;

            var bytes_per_line: usize = 0;
            var bytes_per_line_expected: ?usize = null;

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
                        reader.streamUntilDelimiter(stage.writer(), '\n', options.reserve_len) catch |err| switch (err) {
                            error.EndOfStream => return false,
                            error.StreamTooLong => return error.BufferOverflow,
                            error.NoSpaceLeft => unreachable,
                            else => |e| return e,
                        };
                        try consumer.writeName(std.mem.trim(u8, stage.getWritten(), &std.ascii.whitespace));

                        stage.reset();
                        self.state = .seq;
                    },
                    .seq => blk: {
                        while (true) {
                            switch (reader.readByte() catch |err| switch (err) {
                                error.EndOfStream => return false,
                            }) {
                                '+' => {
                                    self.state = .qual;
                                    try reader.skipUntilDelimiterOrEof('\n');
                                    bases_per_line = 0;
                                    bytes_per_line = 0;
                                    break :blk;
                                },
                                '>' => {
                                    if (self.format == .fastq)
                                        return error.RecordParseError;
                                    self.state = .name; // return back to the name state
                                    return true;
                                },
                                else => |b| {
                                    if (seq_line_count > 0) {
                                        // we're onto the next line need to check line_counts
                                        if (bytes_per_line_expected) |count| {
                                            if (bytes_per_line != count) return error.RecordParseError;
                                        } else bytes_per_line_expected = bytes_per_line;
                                        if (bases_per_line_expected) |count| {
                                            if (bases_per_line != count) return error.RecordParseError;
                                        } else bases_per_line_expected = bases_per_line;
                                    }

                                    bases_per_line = 0;
                                    bytes_per_line = 0;

                                    const app: [1]u8 = .{b};
                                    _ = stage.write(&app) catch |err| switch (err) {
                                        error.NoSpaceLeft => unreachable,
                                    };
                                },
                            }
                            same_line_chunk: while (true) {
                                seq_line_count += 1;
                                reader.streamUntilDelimiter(stage.writer(), '\n', options.reserve_len) catch |err| switch (err) {
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

                                try consumer.writeSeq(bases_buf);

                                seq_number_bases += bases_buf.len;
                                stage.reset();
                                break :same_line_chunk;
                            }
                        }
                    },
                    .qual => {
                        if (qual_line_count > 0) {
                            if (bytes_per_line_expected) |count| {
                                if (bytes_per_line != count) return error.RecordParseError;
                            } else return error.RecordParseError;
                            if (bases_per_line_expected) |count| {
                                if (bases_per_line != count) return error.RecordParseError;
                            } else return error.RecordParseError;

                            bases_per_line = 0;
                            bytes_per_line = 0;
                        }

                        same_line_chunk: while (true) {
                            qual_line_count += 1;
                            reader.streamUntilDelimiter(stage.writer(), '\n', options.reserve_len) catch |err| switch (err) {
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

                            try consumer.writeQual(qual_buf);
                            qual_number_bases += qual_buf.len;
                            stage.reset();
                            break :same_line_chunk;
                        }

                        if (qual_number_bases == seq_number_bases) {
                            if (seq_line_count != qual_line_count) return error.RecordParseError;
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




pub fn fastaReadIterator(stream: anytype, comptime options: FaReadIterOptions) FaReadIterator(@TypeOf(stream), options) {
    return .{ .stream = stream };
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
    var iter = seq.fa.fastaReadIterator(stream, .{});
    const TestCase = struct { name: []const u8, sequence: []const u8, quality: []const u8 };

    // zig fmt: off
    const test_cases =
        [_]TestCase{ 
        .{ 
            .name = "HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG", 
            .sequence = "TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGTATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT", 
            .quality = "HIHIIIIIHIIHGHHIHHIIIIIIIIIIIIIIIHHIIIIIHHIHIIIIIGIHIIIIHHHHHHGHIHIIIIIIIIIIIGHIIIIIGHIIIIHIIHIHHIIIIHIHHIIIIIIIGIIIIIIIHIIIIIGHIIIIHIIIH?DGHEEGHIIIIIIIIIIIHIIHIIIHHIIHIHHIHCHHIIHGIHHHHHHH<GG?B@EHDE-BEHHHII5B@GHHF?CGEHHHDHIHIIH" }, 
        .{ 
            .name = "HWI-D00523:240:HF3WGBCXX:1:1101:5586:3020 1:N:0:CTGTAG", 
            .sequence = "TGGGGAATATTGGGCAATGGGCGGAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCCCTCGGGTTGTAAACTTCTTTTCTATAGGACGAAGAAGTGACGGTACTATAGGAATAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCGAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGAGAGCAAGTCAGATGTGA", 
            .quality = "EHEHGIIIHIGGHGHFEHHEHGCHHGGHIIIGHHIFHHGHHIEHIIIIGHHHHIIGIHGHGGHHHHHCHHHICCHHHHH@HHIIIGCEHGHHGHCHHGDGCGCCEHEEHGIIGHHGHHHIGGFCFHHIHIGIIIHGGHFHIIFEFHIIHIGDHCHFHHGCHCE?GHIIH<C?GHHHHIGFDEHHHHEC88@<@@EHHHIHDH-@HHCDHHDDEHH6@F6@6@@EH@@" }, 
        .{ 
            .name = "HWI-D00523:240:HF3WGBCXX:1:1101:2860:2149 1:N:0:CTGTAG", 
            .sequence = "TAGGGAATATTGCTCAATGGGGGAAACCCTGAAGCAGCAACGCCGCGTGGAGGATGAAGGTTTTAGGATTGTAAACTCCTTTTGTGAGAGAAGATTATGACGGTATCTCACGAATAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCGAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGTAGGCGGGACTGCAAGTTGGGTGTCAAA", 
            .quality = "HGHHGHIIHIIIIIIHHHIIHIIIIIHGHCHIHIIHIIHIIIIIIIIIHIGHIEHIHIIG<FEHHHIHHIIIIIIHIFFHHHHIIIIIHHHHGHFEHHIHHHIHEHHFHHHIHIIIIIIIHIHDHHHHHEHHIIGIIHHIIGHHHIIGDDAGGHHFHHHIHICHHHGH,GHHHHGCEHEG?@6@G?-@>HHHHHHHDEH<@H-@CDD>:E?@GHEF-@E:@H+@-@@" 
        } 
    };
    // zig fmt: on 

    const FaTestRecord = struct {
        const Self = @This();

        name: std.ArrayListUnmanaged(u8),
        qual: std.ArrayListUnmanaged(u8),
        seq: std.ArrayListUnmanaged(u8),

        pub const StreamError = error{OutOfMemory};
        pub const FaStreamConsumer = seq.consumer.SeqQualConsumer(*Self, StreamError, writeName, writeSeq, writeQual);
        fn writeName(self: *Self, buf: []const u8) StreamError!void {
            try self.name.appendSlice(std.testing.allocator, buf);
        }
        fn writeQual(self: *Self, buf: []const u8) StreamError!void {
            try self.qual.appendSlice(std.testing.allocator, buf);
        }
        fn writeSeq(self: *Self, buf: []const u8) StreamError!void {
            try self.seq.appendSlice(std.testing.allocator, buf);
        }
    };
    var faTestRecord: FaTestRecord = .{ .name = .{}, .seq = .{}, .qual = .{} };
    defer {
        faTestRecord.qual.deinit(std.testing.allocator);
        faTestRecord.seq.deinit(std.testing.allocator);
        faTestRecord.name.deinit(std.testing.allocator);
    }

    var con: FaTestRecord.FaStreamConsumer = .{ .context = &faTestRecord };
    for (test_cases) |case| {
        if (try iter.next(con)) {
            try testing.expectEqualStrings(case.name, faTestRecord.name.items);
            try testing.expectEqualStrings(case.sequence, faTestRecord.seq.items);
            try testing.expectEqualStrings(case.quality, faTestRecord.qual.items);
        } else try testing.expect(false);
        faTestRecord.qual.clearRetainingCapacity();
        faTestRecord.seq.clearRetainingCapacity();
        faTestRecord.name.clearRetainingCapacity();
    }
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
