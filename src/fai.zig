const std = @import("std");
const seq = @import("seq.zig");


//faidx – an index enabling random access to FASTA and FASTQ files 
pub fn Entry(comptime container: seq.SeqContianer) type {
    switch(container) {
        .fastq => struct {
            name: [] const u8,
            length: usize,
            offset: usize,
            line_base: usize,
            line_width: usize,
            qual_offset: usize
        },
        .fasta => struct {
            name: [] const u8,
            length: usize,
            offset: usize,
            line_base: usize,
            line_width: usize,
        },
        else => @compileError("Unsupported format: " ++ seq.SeqContianer.name(container))
    }
}

pub const FaiReaderOption =  struct {
    container: seq.SeqContianer, // this is requried
    bufferReserver: usize 
};

pub fn FaiReader(comptime Stream: type, comptime option: FaiReaderOption) type{
   return struct {
        stream: Stream,
        stage: std.BoundedArray(u8, option.bufferReserver),  
        pub const Container = seq.fai.FaiIndex(option.container);
        const Self = @This();

        pub fn next(self: *Self) Container {
            var reader = self.stream.reader();
            var col: usize = 0;
            while(true) {
                const byte: u8 = reader.readByte() catch |err| switch (err) {
                    error.EndOfStream => return null, // we just have a name so the buffer is malformed
                    else => |e| return e,
                };
                switch(byte) {
                    // skip \r
                    std.ascii.control_code.cr => break, 
                    // horizontal tab
                    std.ascii.control_code.ht => {
                        switch(option.container) {
                           .fastq => {
                                switch(col) {
                                    0 => {
                                    }
                                }
                            },
                            .fasta => {

                            },
                            else => unreachable
                        }
                        col += 1;       
                    },
                    // line feed \n
                    std.ascii.control_code.lf => {

                    }
                }
                self.stage.append(byte);
            }

        } 
    };
}


pub fn faiReader(stream: anytype, comptime option: FaiReaderOption) FaiReader(@TypeOf(stream), option){
    return .{}; 
}

test "test read fa records" {
    const file = @embedFile("./test/ce.fa.fai");
    _ = file;

}


