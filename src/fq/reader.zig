const seq = @import("../seq.zig");
const std = @import("std");
const testing = std.testing;

pub const FqReadOptions = struct {};

// read the next record in serial
pub fn readNextFqRecord(reader: anytype, record: anytype) (@TypeOf(reader.*).Error || @TypeOf(record.*).Error || error{ MismatchedQualAndSeq, MalformedIdentifier })!bool {
    record.reset();
    var name_writer = record.nameWriter();
    var seq_writer = record.seqWriter();
    var qual_writer = record.qualWriter();

    var count_seq: u32 = 0;
    var count_qul: u32 = 0;
    var state: enum { start, name, seq, seq_start, qual } = .start;
    blk: while (true) {
        switch (state) {
            .start => {
                const byte: u8 = reader.readByte() catch |err| switch (err) {
                    error.EndOfStream => return false, // we just have a name so the buffer is malformed
                    else => |e| return e,
                };
                state = switch (byte) {
                    '@' => .name,
                    else => return error.MalformedIdentifier,
                };
            },
            .name => {
                const byte: u8 = reader.readByte() catch |err| switch (err) {
                    error.EndOfStream => return false, // we just have a name so the buffer is malformed
                    else => |e| return e,
                };
                state = switch (byte) {
                    '\r' => state, // ignore \r characters
                    '\n' => .seq_start, // end of the line
                    else => l: {
                        try name_writer.writeByte(byte);
                        break :l state;
                    },
                };
            },
            .seq_start => {
                const byte: u8 = reader.readByte() catch |err| switch (err) {
                    error.EndOfStream => return true, // we just have a name so the buffer is malformed
                    else => |e| return e,
                };
                state = switch (byte) {
                    '+' => .qual,
                    else => l: {
                        count_seq += 1;
                        try seq_writer.writeByte(byte);
                        break :l .seq;
                    },
                };
            },
            .seq => {
                const byte: u8 = reader.readByte() catch |err| switch (err) {
                    error.EndOfStream => return false, // we just have a name so the buffer is malformed
                    else => |e| return e,
                };
                state = switch (byte) {
                    '\r' => state,
                    '\n' => .seq_start,
                    else => l: {
                        count_seq += 1;
                        try seq_writer.writeByte(byte);
                        break :l state;
                    },
                };
            },
            .qual => {
                const byte: u8 = reader.readByte() catch |err| switch (err) {
                    error.EndOfStream => return false, // we just have a name so the buffer is malformed
                    else => |e| return e,
                };
                state = switch (byte) {
                    '\r' => state,
                    '\n' => l: {
                        // '@' and '+' can show up in qualifiers
                        //   the number of sequence characters has to equal qualifiers
                        if (count_qul == count_seq)
                            break :blk;
                        if (count_qul > count_seq)
                            return error.MismatchedQualAndSeq;
                        break :l state;
                    },
                    else => l: {
                        count_qul += 1;
                        try qual_writer.writeByte(byte);
                        break :l state;
                    },
                };
            },
        }
    }
    return true;
}


pub fn FqReadIterator(comptime Reader: type, comptime options: FqReadOptions) type {
    _ = options;
    const FqResult = struct {
        seq_offset: usize, // offset in the reader
    };

    return struct {
        reader: Reader,
        offset: u32 = 0,
        const Self = @This();

        pub fn next(self: *Self, record: anytype) (Reader.Error || @TypeOf(record.*).Error || error{ MismatchedQualAndSeq, MalformedIdentifier })! ?FqResult{
            record.reset();
            var name_writer = record.nameWriter();
            var seq_writer = record.seqWriter();
            var qual_writer = record.qualWriter();

            var seq_read_offset: usize = 0;

            var count_seq: u32 = 0;
            var count_qul: u32 = 0;
            var state: enum { start, name, seq, seq_start, qual } = .start;
            finished: while (true) {
                switch (state) {
                    .start => {
                        const byte: u8 = self.reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => return null, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        self.offset += 1;
                        state = switch (byte) {
                            '@' => .name,
                            else => return error.MalformedIdentifier,
                        };
                    },
                    .name => {
                        const byte: u8 = self.reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => return null, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        self.offset += 1;
                        state = switch (byte) {
                            '\r' => state, // ignore \r characters
                            '\n' => .seq_start, // end of the line
                            else => l: {
                                try name_writer.writeByte(byte);
                                break :l state;
                            },
                        };
                    },
                    .seq_start => {
                        const byte: u8 = self.reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => continue :finished, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        self.offset += 1;
                        state = switch (byte) {
                            '+' => .qual,
                            else => l: { 
                                seq_read_offset = self.offset;

                                count_seq += 1;
                                try seq_writer.writeByte(byte);
                                break :l .seq;
                            },
                        };
                    },
                    .seq => {
                        const byte: u8 = self.reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => return null, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        self.offset += 1;
                        state = switch (byte) {
                            '\r' => state,
                            '\n' => .seq_start,
                            else => l: {
                                count_seq += 1;
                                try seq_writer.writeByte(byte);
                                break :l state;
                            },
                        };
                    },
                    .qual => {
                        const byte: u8 = self.reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => return null, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        self.offset += 1;
                        state = switch (byte) {
                            '\r' => state,
                            '\n' => l: {
                                // '@' and '+' can show up in qualifiers
                                //   the number of sequence characters has to equal qualifiers
                                if (count_qul == count_seq)
                                    break :finished;
                                if (count_qul > count_seq)
                                    return error.MismatchedQualAndSeq;
                                break :l state;
                            },
                            else => l: {
                                count_qul += 1;
                                try qual_writer.writeByte(byte);
                                break :l state;
                            },
                        };
                    },
                }
            }
            return .{
                .seq_offset = seq_read_offset
            };
        }
    };
}

test "read fq records serially" {
    const file = @embedFile("./test/t1.fq");
    var stream = std.io.fixedBufferStream(file);

    var current = seq.fq.FqArrayList.init(std.testing.allocator);
    var fqReadIter = seq.fq.createReader(stream.reader(), .{});
    // var iter = fq_slice_iter(file);
    if(try fqReadIter.next(&current)) |_| {
        try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG", current.name.items);
        try testing.expectEqualStrings("TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGTATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT", current.seq.items);
    }

//    try testing.expectEqual(try readNextFqRecord(&reader, &current), true);
//    try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:5586:3020 1:N:0:CTGTAG", current.name.items);
//
//    try testing.expectEqual(try readNextFqRecord(&reader, &current), true);
//    try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:2860:2149 1:N:0:CTGTAG", current.name.items);

    current.deinit();
}
