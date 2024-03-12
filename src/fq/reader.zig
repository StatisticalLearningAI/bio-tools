const seq = @import("../seq.zig");
const std = @import("std");
const testing = std.testing;

pub const Option = struct {
    skipQual: bool = false, // skip qualifiers
};


pub fn FqReader(comptime Stream: type, comptime options: Option) type {

    const FqResult = struct {
        seq_offset: usize, // offset in the reader
        sequence_len: usize, // the number of sequence length
        slice: seq.fq.FqSlice
    };

    return struct {
        stream: Stream,
        offset: u32 = 0,
        chunk: std.ArrayList(u8),  

        const Self = @This();

        pub fn deinit(self: *Self) void {
            self.chunk.deinit();
        }
        // move seek to the next the record
        pub fn seekAndReadUsingFai(self: *Self, lookup: seq.fai.FaiSlice) (Stream.SeekableStream.Error)!?FqResult{
            _ = lookup;
            var seek = self.stream.seekableStream();
            _ = seek;
            
            return .{

            };
        }

        pub fn readNext(self: *Self) (Stream.Reader.Error || error{ OutOfMemory, MismatchedQualAndSeq, MalformedIdentifier })!?FqResult {
            self.chunk.shrinkRetainingCapacity(0);
            var seq_read_offset: usize = 0;

            var reader = self.stream.reader();

            var count_seq: u32 = 0;
            var count_qul: u32 = 0;
            var state: enum { start, name, seq, seq_start, qual } = .start;
            
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
                        self.offset += 1;
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
                        self.offset += 1;
                        state = switch (byte) {
                            '\r' => state, // ignore \r characters
                            '\n' => .seq_start, // end of the line
                            else => l: {
                                try self.chunk.append(byte);
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
                        self.offset += 1;
                        state = switch (byte) {
                            '+' => .qual,
                            else => l: {
                                seq_read_offset = self.offset;

                                count_seq += 1;
                                try self.chunk.append(byte);
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
                        self.offset += 1;
                        state = switch (byte) {
                            '\r' => state,
                            '\n' => .seq_start,
                            else => l: {
                                count_seq += 1;
                                try self.chunk.append(byte);
                                seq_len += 1;
                                break :l state;
                            },
                        };
                    },
                    .qual => {
                        const byte: u8 = reader.readByte() catch |err| switch (err) {
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
                                //try qual_writer.writeByte(byte);
                                if (!options.skipQual) {
                                    try self.chunk.append(byte);
                                    qual_len += 1; 
                                }
                                break :l state;
                            },
                        };
                    },
                }
            }
            return .{ .seq_offset = seq_read_offset, .sequence_len = count_seq, .slice = .{
                .name = self.chunk.items[0..name_len],
                .seq = self.chunk.items[name_len..name_len + seq_len],
                .qual = self.chunk.items[name_len + seq_len..qual_len + seq_len + name_len]
            } };
        }
    };
}

test "read fq records serially" {
    const file = @embedFile("./test/t1.fq");

    var stream = std.io.fixedBufferStream(file);
    var fqReadIter = seq.fq.createReader(std.testing.allocator,stream, .{});
    defer fqReadIter.deinit();
    if (try fqReadIter.readNext()) |result| {
        try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG", result.slice.name);
        try testing.expectEqualStrings("TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGTATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT", result.slice.seq);
    } else try testing.expect(false);

    if (try fqReadIter.readNext()) |result| {
        try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:5586:3020 1:N:0:CTGTAG", result.slice.name);
    } else try testing.expect(false);

    if (try fqReadIter.readNext()) |result| {
        try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:2860:2149 1:N:0:CTGTAG", result.slice.name);
    } else try testing.expect(false);


}
