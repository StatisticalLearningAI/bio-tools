

pub fn readNextFaiRecord(reader: anytype, record: anytype) void {
    record.reset();

    var name_writer = reader.nameWriter();
    _ = name_writer;

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
