const std = @import("std");
const testing = std.testing;


pub const SeqContianer = @import("container.zig").SeqContianer; 
pub const fai = @import("./fai.zig"); 
pub const fa = @import("./fa.zig"); 

//pub const fqWriter = fa.fqWriter;
//pub const fqReader = fa.fqReader;

//pub const faiWriter = fai.faiWriter;
pub const faiReader = fai.faiReaderIterator;

test {
    _ = fai;
    _ = fa;
}

