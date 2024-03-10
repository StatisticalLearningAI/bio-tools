const std = @import("std");
const Block = struct { begin: u32, end: u32, next: u32 };

fn str_next_line(index: usize, slice: []const u8) usize {
    var idx = index;
    while (idx < slice.len) : (idx += 1) {
        if (slice[idx] == '\n') {
            idx += 1;
            return idx;
        } else if (slice[idx] == '\r') {
            idx += 1;
            if (idx < slice.len and slice[idx] == '\n') {
                idx += 1;
            }
            return idx;
        }
    }
    return slice.len;
}

fn str_r_trim(slice: []const u8) []const u8 {
    var result = slice;
    while (result.len > 0) {
        if (result[result.len - 1] == '\r' or result[result.len - 1] == '\n') {
            result = result[0..(result.len - 1)];
        }
        break;
    }
    return result;
}

inline fn is_carriage_return_char(c: u8) bool {
    return c == '\r' or c == '\n';
}

inline fn is_empty_char(c: u8) bool {
    return c == '\r' or c == '\n' or c == ' ';
}

fn stream_until_carriage_return(reader: anytype, writer: anytype) !void {
    while (true) {
        const byte: u8 = try reader.readByte();
        if (is_carriage_return_char(byte)) return;
        try writer.writeByte(byte);
    }
}

fn str_trim(slice: []const u8) []const u8 {
    if (slice.len == 0) return slice;
    var left: usize = 0;
    var right: usize = slice.len - 1;
    while (left < right) {
        if (is_empty_char(slice[left])) {
            left += 1;
            continue;
        }
        if (is_empty_char(slice[right])) {
            right -= 1;
            continue;
        }
        break;
    }
    return slice[left .. right + 1];
}

fn str_start_with_ignore_spaces(slice: []const u8, c: u8) ?u32 {
    var index: u32 = 0;
    while (index < slice.len) : (index += 1) {
        if (slice[index] == ' ') continue;
        if (slice[index] == c)
            return index;
        break;
    }
    return null;
}

const FqRowType = enum { None, Seq, Quality };
fn fq_line_type(slice: []const u8, offset: ?*usize) FqRowType {
    if (str_start_with_ignore_spaces(slice, '@')) |o| {
        if (offset) |of| (of.*) += o;
        return FqRowType.Seq;
    } else if (str_start_with_ignore_spaces(slice, '+')) |o| {
        if (offset) |of| (of.*) += o;
        return FqRowType.Quality;
    }
    return FqRowType.None;
}

fn to_record(name_: ?[]const u8, seq_: ?[]const u8, qual_: ?[]const u8) ?FqRecord {
    if (seq_) |s| {
        return .{ .name = name_, .seq = s, .qual = qual_ };
    }
    return null;
}


pub const FqRecord = struct {
    pub const FqArray = std.ArrayListAligned(u8, 16);

    name: FqArray,
    qual: FqArray,
    seq: FqArray,

    pub fn clear(self: *FqRecord) void {
        self.name.resize(0) catch {};
        self.qual.resize(0) catch {};
        self.seq.resize(0) catch {};
    }

    // qual and sequences should match 
    pub fn isQualAndSeqValid(self: FqRecord) bool {
        return self.qual.len > 0 and self.qual.len == self.seq.len;
    }

    pub fn deinit(self: *FqRecord) void {
        self.name.deinit();
        self.qual.deinit();
        self.seq.deinit();
    }

    pub fn init(allocator: std.mem.Allocator) FqRecord {
        return .{
            .name = FqArray.init(allocator),
            .qual = FqArray.init(allocator),
            .seq = FqArray.init(allocator),
        };
    }
};


pub fn FqWriter(comptime Writer: type) type {
    _ = Writer;
    return struct {

    };
}

pub fn FqReader(comptime Reader: type) type {
    return struct {
        reader: Reader,
        start_read_name: bool,

        const Self = @This();

        pub fn read_next_record(self: *Self, rec: *FqRecord) (Reader.Error || error{ OutOfMemory, EndOfRecord })!void {
            rec.clear();
            var qual = false;
            while (self.start_read_name) {
                const b2: u8 = self.reader.readByte() catch |err| switch (err) {
                    error.EndOfStream => return error.EndOfRecord,
                    else => |e| return e,
                };
                if (is_carriage_return_char(b2)) break;
                rec.name.append(b2) catch |err| return err;
            }
            self.start_read_name = false;

            while (true) {
                const byte: u8 = self.reader.readByte() catch |err| switch (err) {
                    error.EndOfStream => return error.EndOfRecord,
                    else => |e| return e,
                };
                if (is_empty_char(byte)) continue;
                switch (byte) {
                    '@' => {
                        if (rec.seq.items.len > 0) {
                            self.start_read_name = true;
                            return; // we found a new sequence so terminate early
                        }
                        while (true) {
                            const b2: u8 = self.reader.readByte() catch |err| switch (err) {
                                error.EndOfStream => return,
                                else => |e| return e,
                            };
                            if (is_carriage_return_char(b2)) break;
                            rec.name.append(b2) catch |err| return err;
                        }
                    },
                    '+' => {
                        qual = true;
                        try self.reader.skipUntilDelimiterOrEof('\n');
                    },
                    else => {},
                }
                while (true) {
                    const b2 = self.reader.readByte() catch |err| switch (err) {
                        error.EndOfStream => return,
                        else => |e| return e,
                    };
                    if (is_carriage_return_char(b2)) break;
                    if (qual) {
                        rec.qual.append(b2) catch |err| return err;
                    } else {
                        rec.seq.append(b2) catch |err| return err;
                    }
                }
            }
        }
    };
}

pub fn fq_stream_reader(reader: anytype) FqReader(@TypeOf(reader)) {
    return .{ .reader = reader, .start_read_name = false};
}

const testing = std.testing;
test "fq_read_test" {
    const file = @embedFile("./files/t1.fq");
    var stream = std.io.fixedBufferStream(file);

    var current = FqRecord.init(std.testing.allocator);
    var fq = fq_stream_reader(stream.reader());
    // var iter = fq_slice_iter(file);
    try fq.read_next_record(&current);
    try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG", current.name.items);
    try fq.read_next_record(&current);
    try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:5586:3020 1:N:0:CTGTAG", current.name.items);

    current.deinit();
}
