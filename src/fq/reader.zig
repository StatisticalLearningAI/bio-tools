const seq = @import("../seq.zig"); 
const std = @import("std");
const testing = std.testing;

pub fn FqReader(comptime Reader: type) type {
    return struct {
        reader: Reader,
        start_read_name: bool,
        const Self = @This();
    
        pub fn read_next_record(self: *Self, writer: anytype) (Reader.Error || error{ OutOfMemory, EndOfRecord })!void {
            writer.resetRecord();
            var qual = false;
            var markStart = false;
            while (self.start_read_name) {
                markStart = true;
                const b2: u8 = self.reader.readByte() catch |err| switch (err) {
                    error.EndOfStream => return error.EndOfRecord,
                    else => |e| return e,
                };
                if (seq.str.is_carriage_return_char(b2)) break;
                writer.name.writeByte(b2) catch |err| return err;
            }
            self.start_read_name = false;

            while (true) {
                const byte: u8 = self.reader.readByte() catch |err| switch (err) {
                    error.EndOfStream => return error.EndOfRecord,
                    else => |e| return e,
                };
                if (seq.str.is_empty_char(byte)) continue;
                switch (byte) {
                    '@' => {
                        if (markStart == true) {
                            self.start_read_name = true;
                            return; // we found a new sequence so terminate early
                        }
                        markStart = true;
                        while (true) {
                            const b2: u8 = self.reader.readByte() catch |err| switch (err) {
                                error.EndOfStream => return,
                                else => |e| return e,
                            };
                            if (seq.str.is_carriage_return_char(b2)) break;
                            writer.name.writeByte(b2) catch |err| return err;
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
                    if (seq.str.is_carriage_return_char(b2)) break;
                    if (qual) {
                        writer.qual.writeByte(b2) catch |err| return err;
                    } else {
                        writer.seq.writeByte(b2) catch |err| return err;
                    }
                }
            }
        }
    };
}

test "fq_read_test" {
    const file = @embedFile("./test/t1.fq");
    var stream = std.io.fixedBufferStream(file);

    var current = seq.fq.FqArrayListRecord.init(std.testing.allocator);
    var writer = seq.fq.createReader(stream.reader());
    // var iter = fq_slice_iter(file);
    try writer.read_next_record(&current.writer());
    try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG", current.name.items);
    try writer.read_next_record(&current.writer());
    try testing.expectEqualStrings("HWI-D00523:240:HF3WGBCXX:1:1101:5586:3020 1:N:0:CTGTAG", current.name.items);
    try writer.read_next_record(&current.writer());
    try testing.expectEqualStrings("HHIIIGCEHGHHGHCHHGDGCGCCEHEEHGIIGHHGHHHIGGFCFHHIHIGIIIHGGHFHIIFEFHIIHIGDHCHFHH", current.name.items);

    current.deinit();
}

