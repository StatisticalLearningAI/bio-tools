const std = @import("std");


const Record = struct {

};

pub fn load_buf(slice: []const u8) !void {
   var cur = 0;
   while(cur < slice.len) {
        if(slice[cur] == '>') {

            const Test = struct {
                v: u8,
                k: u32
            };
            const vv: Test = &slice[cur];
            vv.v = 32;

            
            cur += 1;
            while(cur < slice.len) {

            }
        }
    }
}
