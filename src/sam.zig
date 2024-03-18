const seq = @import("./seq.zig");
const std = @import("std");
const testing = std.testing;

// https://samtools.github.io/hts-specs/SAMv1.pdf

pub const SortOrder = enum {
    Unknown, 
    UnSorted,
    QueryName,
    Coordinate
};

pub const FieldType = enum { query_name, flags, reference_name, pos };

pub const HDTagType = enum { version, sorted_alignment };

pub const AlignmentColumns = {};

const HDTag = union(HDTagType) {
    version: struct { major: u16, minor: u16 },
};


pub const ReadGroup = struct {};

pub const RefSeqName = struct {};

pub const SamAlignmentFlags = struct {
    paired: bool = false, // 0x1 the read is paired in sequencing, no matter whether it is mapped in a pair
    proper_pair: bool = false, // 0x2 each segment properly aligned according to the aligner
    unmap: bool = false, // 0x4 segment is unmapped conflicts with proper_pair
    munmap: bool = false, // 0x8 next segment in the template unmapped
    reverse: bool = false, // 0x10 SEQ being reverse complemented
    mate_reverse: bool = false, // 0x20 SEQ of the next segment in the template being reverse complemented
    read1: bool = false, // 0x40 first segment in the template
    read2: bool = false, // 0x80 second segment in the template
    secondary: bool = false, // 0x100 secondary alignment
    qc_failure: bool = false, // 0x200 not passing filters, such as platform/vendor quality controls QC
    dup: bool = false, // 0x400 PCR or optical duplicate
    sup_alignment: bool = false, // 0x800 supplementary allignment

    pub const PariedBit = 0x1;
    pub const ProperPairedBit = 0x2;
    pub const Unmapped = 0x4;
    pub const MateUnmapped = 0x8;
    pub const Reverse = 0x10;
    pub const MateReverse = 0x20;
    pub const Read1 = 0x40;
    pub const Read2 = 0x80;
    pub const Secondary = 0x100;
    pub const QCFailure = 0x200;
    pub const Dup = 0x400;
    pub const SupplementaryAlignment = 0x800;
    
    pub fn toSamAlignmentFlags(value: anytype) SamAlignmentFlags {
        return .{
            .paired = (SamAlignmentFlags.PariedBit & value) > 0,  
            .proper_pair = (SamAlignmentFlags.ProperPairedBit & value) > 0,
            .unmap = (SamAlignmentFlags.Unmapped & value) > 0,
            .munmap = (SamAlignmentFlags.MateUnmapped & value) > 0,
            .reverse = (SamAlignmentFlags.Reverse & value) > 0,
            .mate_reverse = (SamAlignmentFlags.MateReverse & value) > 0,
            .read1 = (SamAlignmentFlags.Read1 & value) > 0,
            .read2 = (SamAlignmentFlags.Read2 & value) > 0,
            .secondary = (SamAlignmentFlags.Secondary & value) > 0,
            .qc_failure = (SamAlignmentFlags.QCFailure & value) > 0,
            .dup = (SamAlignmentFlags.Dup & value) > 0,
            .sup_alignment = (SamAlignmentFlags.SupplementaryAlignment & value) > 0
        };
    }

};

pub const ReaderOptions = struct {
    reserve_len: usize = 4096,
};

pub fn SamReadIterator(comptime Stream: type, comptime options: ReaderOptions) type {
    return struct {
        const Self = @This();
        stream: Stream,

        state: enum {
            start,
            begin_alignment,
            alignment
        } = .start,
        buffer: [options.reserve_len]u8 = undefined,
        allocator: std.mem.Allocator,
        pos: usize = 0,

        // collected from the header @HD
        version: struct { major: u32, minor: u32 },
        sort_order: SortOrder, 

        // read group @RG
        read_groups: std.ArrayListUnmanaged(ReadGroup),
        // list of refences @SQ
        reference_seqs: std.ArrayListUnmanaged(RefSeqName),

        // init and read the header porition of the sam
        pub fn init(self: *Self, allocator: std.mem.Allocator) void {
            self.allocator = allocator;
        }

        pub fn deinit(self: *Self) void {
            _ = self;
        }

        pub fn parseHeader(self: *Self) void {
            switch(self.state) {
                .start => {
                    var reader = self.stream.reader();
                    const byte: u8 = reader.readByte() catch |err| switch (err) {
                        error.EndOfStream => return null, // we just have a name so the buffer is malformed
                        else => |e| return e,
                    };
                    switch(byte) {
                        // horizontal tab
                        std.ascii.control_code.ht, std.ascii.control_code.lf => l: {
                            if(self.buffer[0] == '@') {

                                if(std.mem.eql(u8, self.buffer[0..self.pos], "@HD")) {
                                }
                            } else {
                                self.state = .begin_alignment;
                                break :l;
                            }
                        },
                        else => {
                            self.buffer[self.pos] = byte;
                            self.pos += 1;
                        }
                    }
                }
            }

        }

        pub fn next(self: *Self, sequence: anytype) struct {
            name: []const u8,
            flags: struct {},
        } {
            std.debug.assert(self.state == .alignment); // the header has to be parsed first
            


            _ = sequence;
        }
    };
}

//1 0x1 template having multiple segments in sequencing
//2 0x2 each segment properly aligned according to the aligner
//4 0x4 segment unmapped
//8 0x8 next segment in the template unmapped
//16 0x10 SEQ being reverse complemented
//32 0x20 SEQ of the next segment in the template being reverse complemented
//64 0x40 the first segment in the template
//128 0x80 the last segment in the template
//256 0x100 secondary alignment
//512 0x200 not passing filters, such as platform/vendor quality controls
//1024 0x400 PCR or optical duplicate
//2048 0x800 supplementary alignment

//1 QNAME String [!-?A-~]{1,254} Query template NAME
//2 FLAG Int [0, 216 − 1] bitwise FLAG
//3 RNAME String \*|[:rname:∧*=][:rname:]* Reference sequence NAME11
//4 POS Int [0, 231 − 1] 1-based leftmost mapping POSition
//5 MAPQ Int [0, 28 − 1] MAPping Quality
//6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
//7 RNEXT String \*|=|[:rname:∧*=][:rname:]* Reference name of the mate/next read
//8 PNEXT Int [0, 231 − 1] Position of the mate/next read
//9 TLEN Int [−231 + 1, 231 − 1] observed Template LENgth
//10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
//11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33

test "read sam record" {
    const file = @embedFile("./test/t1.fq");
    var stream = std.io.fixedBufferStream(file);
    var iter = seq.fa.fastaReadIterator(stream, .{});
    _ = iter;

    const TestResult = struct {
        const Self = @This();
        name: std.ArrayListUnmanaged(u8),
        qual: std.ArrayListUnmanaged(u8),
        seq: std.ArrayListUnmanaged(u8),
    };
    var faTestRecord: TestResult = .{ .name = .{}, .seq = .{}, .qual = .{} };

    var consumer = struct {
        context: *TestResult,
        const Self = @This();
        pub const Error = error{OutOfMemory};
        pub inline fn writeSeq(buf: []const u8) void {
            _ = buf;
        }
        pub inline fn writeQual(buf: []const u8) void {
            _ = buf;
        }
    }{ .context = &faTestRecord };
    _ = consumer;
}
