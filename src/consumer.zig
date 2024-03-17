pub fn SeqConsumer (comptime Context: type, 
        comptime StreamError: type,
        comptime fnWriteName: fn (self: Context, buf: []const u8) StreamError!void, 
        comptime fnWriteSeq: fn (self: Context, buf: []const u8) StreamError!void) type {
    return struct {
        context: Context,
        pub const Self = @This();
        pub const Error = StreamError;
        pub inline fn writeName(self: Self, buf: []const u8) StreamError!void {
            try fnWriteName(self.context, buf);
        }
        pub inline fn writeSeq(self: Self, buf: []const u8) StreamError!void {
            try fnWriteSeq(self.context, buf);
        }
    };
}

pub fn SeqQualConsumer (comptime Context: type, comptime StreamError: type, comptime fnWriteName: fn (self: Context, buf: []const u8) StreamError!void, comptime fnWriteSeq: fn (self: Context, buf: []const u8) StreamError!void, comptime fnWriteQual: fn (self: Context, buf: []const u8) StreamError!void) type {
    return struct {
        context: Context,
        pub const Self = @This();
        pub const Error = StreamError;
        pub inline fn writeName(self: Self, buf: []const u8) StreamError!void {
            try fnWriteName(self.context, buf);
        }
        pub inline fn writeSeq(self: Self, buf: []const u8) StreamError!void {
            try fnWriteSeq(self.context, buf);
        }
        pub inline fn writeQual(self: Self, buf: []const u8) StreamError!void {
            try fnWriteQual(self.context, buf);
        }
    };
}

