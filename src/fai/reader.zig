const seq = @import("../seq.zig");

pub fn FaiReader(comptime Stream: type) type{
   return struct {
        stream: Stream,
        
        const Self = @This();
    
        pub fn deinit(self: *Self) void {
           _ = self;

        }

        pub fn readNex(self: *Self) seq.fai.FaiSlice {
           _ = self;
            
        } 
    };
}


//pub fn FaiReader(comptime Reader: type) type {
//    return struct {
//        const Self = @This();
//        
//        reader: Reader,
//        
//
//        pub fn readNextRecord(self: *Self) bool{
//
//
//            return true;
//        }
//    };
//}
