const std = @import("std");

const Block = struct { begin: u32, end: u32, next: u32 };

inline fn is_carriage_return(c: u8) bool {
    return c == '\r' or c == '\n';
}

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
fn resolve_line_type(slice: []const u8, offset: ?*usize) FqRowType {
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

const FqResult = enum {
    FqSuccess,
    FqMismatchedQualLength,
};

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

    pub fn validate(self: *FqRecord) FqResult {
        _ = self;
        return FqResult.FqSuccess;
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

const FqSliceIter = struct {
    index: usize,
    slice: []const u8,
    pub fn read_next_record(self: *FqSliceIter, rec: *FqRecord) ?*FqRecord {
        rec.clear();
        var qual = false;
        main: while (self.index < self.slice.len) {
            switch (resolve_line_type(self.slice[self.index..], &self.index)) {
                FqRowType.Quality => {
                    qual = true;
                    self.index = str_next_line(self.index, self.slice);
                },
                FqRowType.Seq => {
                    if (rec.seq.items.len > 0)
                        break :main;
                    const next_idx = str_next_line((self.index + 1), self.slice);
                    rec.name.appendSlice(str_r_trim(self.slice[(self.index + 1)..next_idx])) catch {};
                    self.index = next_idx;
                },
                FqRowType.None => {
                    const next_idx = str_next_line(self.index, self.slice);
                    if (qual) {
                        rec.qual.appendSlice(str_r_trim(self.slice[self.index..next_idx])) catch {};
                    } else {
                        rec.seq.appendSlice(str_r_trim(self.slice[self.index..next_idx])) catch {};
                    }
                    self.index = next_idx;
                },
            }
        }
        return rec;
    }

    pub fn init(slice: []const u8) FqSliceIter {
        return .{ .index = 0, .slice = slice };
    }
};

const testing = std.testing;
test "fq_read_test" {
    const file = @embedFile("./files/t1.fq");
    var current = FqRecord.init(std.testing.allocator);
    var iter = FqSliceIter.init(file);
    if (iter.read_next_record(&current)) |rec| {
        try testing.expectEqualStrings(rec.name.items, "HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG");
    }
    if (iter.read_next_record(&current)) |rec| {
        try testing.expectEqualStrings(rec.name.items, "HWI-D00523:240:HF3WGBCXX:1:1101:5586:3020 1:N:0:CTGTAG");
    }
    current.deinit();
}

