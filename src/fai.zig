const std = @import("std");
const seq = @import("seq.zig");
const testing = std.testing;

const ColumnIdx = enum(u8) { name = 0, length = 1, offset = 2, line_bases = 3, line_width = 4, qual_offset = 5 };

//faidx – an index enabling random access to FASTA and FASTQ files
// if the reserved is 0 then we assume its a const []u8 else the entry takes up some amount of reserved memory
pub fn FaiEntry(comptime container: seq.SeqContianer, comptime reserved_len: usize) type {
    return switch (container) {
        .fastq => struct {
            pub const NumberOfColumns = 6;

            name: if (reserved_len == 0) []const u8 else [reserved_len:0]u8,
            length: u32,
            offset: u64,
            line_base: u32,
            line_width: u32,
            qual_offset: u64,
        },
        .fasta => struct {
            pub const NumberOfColumns = 5;

            name: if (reserved_len == 0) []const u8 else [reserved_len:0]u8,
            length: u32,
            offset: u64,
            line_base: u32,
            line_width: u32,
        },
        else => @compileError("Unsupported format: " ++ seq.SeqContianer.name(container)),
    };
}

pub fn FaiTaggedEntry(comptime reserve_len: usize) type {
    return union(seq.SeqContianer) { fastq: FaiEntry(.fastq, reserve_len), fasta: FaiEntry(.fasta, reserve_len) };
}

pub const FaiReaderOption = struct {
    container: seq.SeqContianer, // this is requried
    bufferReserver: usize = 2048,
};

pub fn FaiReader(comptime Stream: type, comptime option: FaiReaderOption) type {
    return struct {
        stream: Stream,
        // first part for the name rest is reserved for just for parsing values
        // |  name | .. parsing .. |
        stage: std.BoundedArray(u8, option.bufferReserver),

        pub const Iterable = FaiEntry(option.container, 0);

        const Self = @This();
        pub fn next(self: *Self) (error{ UnexpectedNumColumn, InvalidCharacter, Overflow })!?Iterable {
            var reader = self.stream.reader();
            var col: u8 = 0;

            var name_len: usize = 0; // reserved the first chunk of the buffer for the name
            var length: u32 = 0;
            var offset: u64 = 0;
            var line_base: u32 = 0;
            var line_width: u32 = 0;
            var qual_offset: u64 = 0;

            blk: while (true) {
                const byte: u8 = reader.readByte() catch |err| switch (err) {
                    error.EndOfStream => return null, // the buffer is exausted 
                    else => |e| return e,
                };
                switch (byte) {
                    // skip \r
                    std.ascii.control_code.cr => break,
                    // horizontal tab
                    std.ascii.control_code.ht => {
                        switch (comptime option.container) {
                            .fastq => {
                                switch (@as(ColumnIdx, @enumFromInt(col))) {
                                    ColumnIdx.name => name_len = self.stage.len,
                                    ColumnIdx.length => length = try std.fmt.parseUnsigned(u32, self.stage.slice()[name_len..], 10),
                                    ColumnIdx.offset => offset = try std.fmt.parseUnsigned(u64, self.stage.slice()[name_len..], 10),
                                    ColumnIdx.line_bases => line_base = try std.fmt.parseUnsigned(u32, self.stage.slice()[name_len..], 10),
                                    ColumnIdx.line_width => line_width = try std.fmt.parseUnsigned(u32, self.stage.slice()[name_len..], 10),
                                    ColumnIdx.qual_offset => qual_offset = try std.fmt.parseUnsigned(u64, self.stage.slice()[name_len..], 10),
                                    else => unreachable
                                }
                            },
                            .fasta => {
                                switch (@as(ColumnIdx, @enumFromInt(col))) {
                                    ColumnIdx.name => name_len = self.stage.len,
                                    ColumnIdx.length => length = try std.fmt.parseUnsigned(u32, self.stage.slice()[name_len..], 10),
                                    ColumnIdx.offset => offset = try std.fmt.parseUnsigned(u64, self.stage.slice()[name_len..], 10),
                                    ColumnIdx.line_bases => line_base = try std.fmt.parseUnsigned(u32, self.stage.slice()[name_len..], 10),
                                    ColumnIdx.line_width => line_width = try std.fmt.parseUnsigned(u32, self.stage.slice()[name_len..], 10),
                                    else => unreachable
                                }
                            },
                            else => unreachable,
                        }
                        try self.stage.resize(name_len);
                        col += 1;
                    },
                    // line feed \n
                    std.ascii.control_code.lf => {
                        break :blk;
                    },
                    else => try self.stage.append(byte),
                }
            }
            
            return switch (comptime option.container) {
                .fastq => .{
                    .name = self.stage.slice()[0..name_len],
                    .length = length,
                    .offset = offset,
                    .line_base = line_base,
                    .line_width = line_width,
                    .qual_offset = qual_offset 
                },
                .fasta => .{
                    .name = self.stage.slice()[0..name_len],
                    .length = length,
                    .offset = offset,
                    .line_base = line_base,
                    .line_width = line_width,
                },
                else => unreachable
            };
        }
    };
}

pub fn faiReader(stream: anytype, comptime option: FaiReaderOption) FaiReader(@TypeOf(stream), option) {
    return .{ .stream = stream, .stage = .{} };
}

test "test read fa records" {
    const file = @embedFile("./test/ce.fa.fai");
    var stream = std.io.fixedBufferStream(file);
    var reader = faiReader(stream, .{ .container = seq.SeqContianer.fasta });

    if (try reader.next()) |result| {
        try testing.expectEqual(@as(u64, 1009800), result.length);
    } else try testing.expect(false);
}
