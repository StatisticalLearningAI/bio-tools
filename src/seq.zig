const std = @import("std");
const testing = std.testing;


pub const SeqContianer = @import("container.zig").SeqContianer; 
pub const fq = @import("./fq.zig"); 
pub const fai = @import("./fai.zig"); 
pub const fa = @import("./fa.zig"); 

pub const fqWriter = fq.fqWriter;
pub const fqReader = fq.fqReader;

//pub const faiWriter = fai.faiWriter;
pub const faiReader = fai.faiReader;

test {
    _ = fq;
    _ = fai;
}

