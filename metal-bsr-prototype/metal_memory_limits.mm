#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <cstdio>

// Query allocation limits without creating any GPU buffers or dispatching work.
int main(){@autoreleasepool{
    id<MTLDevice> device=MTLCreateSystemDefaultDevice();
    if(!device){std::fputs("No Metal device\n",stderr);return 1;}
    NSDictionary* values=@{
        @"device":device.name,
        @"physical_memory_bytes":@([NSProcessInfo processInfo].physicalMemory),
        @"has_unified_memory":@(device.hasUnifiedMemory),
        @"recommended_max_working_set_bytes":@(device.recommendedMaxWorkingSetSize),
        @"max_buffer_length_bytes":@(device.maxBufferLength)
    };
    NSError* error=nil;
    NSData* json=[NSJSONSerialization dataWithJSONObject:values options:NSJSONWritingPrettyPrinted error:&error];
    if(!json){std::fprintf(stderr,"JSON: %s\n",error.localizedDescription.UTF8String);return 1;}
    std::fwrite(json.bytes,1,json.length,stdout);std::fputc('\n',stdout);return 0;
}}
