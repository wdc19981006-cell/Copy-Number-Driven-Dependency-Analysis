"""Record native failure location in an owned diagnostic R process."""
from pathlib import Path
import ctypes as c,subprocess,struct,json
ROOT=Path(__file__).resolve().parents[2]
k=c.WinDLL('kernel32',use_last_error=True);ps=c.WinDLL('psapi')
class Event(c.Structure):_fields_=[('code',c.c_uint32),('pid',c.c_uint32),('tid',c.c_uint32),('data',c.c_uint64*20)]
k.WaitForDebugEvent.argtypes=[c.POINTER(Event),c.c_uint32];k.ContinueDebugEvent.argtypes=[c.c_uint32,c.c_uint32,c.c_uint32]
ps.GetMappedFileNameW.argtypes=[c.c_void_p,c.c_void_p,c.c_wchar_p,c.c_uint32]
process=subprocess.Popen(['D:/R/R-4.5.0/bin/Rscript.exe','--vanilla','-e','library(cli);cat("Done\\n")'],cwd=ROOT,creationflags=1,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
handle=None;modules=[];fail=[];live=set();handles={};module_map={};exits={};threads={}
k.GetThreadContext.argtypes=[c.c_void_p,c.c_void_p];k.ReadProcessMemory.argtypes=[c.c_void_p,c.c_void_p,c.c_void_p,c.c_size_t,c.c_void_p]
while True:
 e=Event()
 if not k.WaitForDebugEvent(c.byref(e),30000):break
 raw=bytes(e.data)
 if e.code==3:
  live.add(e.pid);handles[e.pid]=struct.unpack_from('Q',raw,8)[0];module_map[e.pid]=[]
  threads[e.tid]=struct.unpack_from('Q',raw,16)[0]
 if e.code==2:threads[e.tid]=struct.unpack_from('Q',raw,0)[0]
 handle=handles.get(e.pid)
 if e.code==6:
  base=struct.unpack_from('Q',raw,8)[0];name=c.create_unicode_buffer(4096)
  if handle:ps.GetMappedFileNameW(handle,base,name,4096)
  module_map[e.pid].append({'base':base,'path':name.value})
 status=0x00010002
 if e.code==1:
  exception=struct.unpack_from('I',raw,0)[0]
  if exception in (0xC0000005,0xC000001D,0xC0000409):
   address=struct.unpack_from('Q',raw,16)[0];mod=max((m for m in module_map[e.pid] if m['base']<=address),key=lambda m:m['base'],default={})
   fail.append({'exception':hex(exception),'address':hex(address),'module':mod.get('path'),'offset':hex(address-mod.get('base',address)),'first_chance':struct.unpack_from('I',raw,152)[0]})
   fail[-1]['memory_target']=hex(struct.unpack_from('Q',raw,40)[0])
   context=c.create_string_buffer(1248);pointer=c.addressof(context);pointer=(pointer+15)&~15
   c.c_uint32.from_address(pointer+48).value=0x100003
   if k.GetThreadContext(threads[e.tid],pointer):
    rsp=c.c_uint64.from_address(pointer+152).value;stack=c.create_string_buffer(1024);received=c.c_size_t()
    k.ReadProcessMemory(handle,rsp,stack,1024,c.byref(received));frames=[]
    for value in struct.unpack('<128Q',stack.raw):
     m=max((m for m in module_map[e.pid] if m['base']<=value),key=lambda m:m['base'],default={});offset=value-m.get('base',value)
     if m and offset<0x3000000:frames.append({'module':m['path'],'offset':hex(offset)})
    fail[-1]['stack_candidates']=frames[:35]
  if exception!=0x80000003:status=0x80010001
 k.ContinueDebugEvent(e.pid,e.tid,status)
 if e.code==5:
  exits[e.pid]=struct.unpack_from('I',raw,0)[0];live.discard(e.pid)
  if not live:break
out,err=process.communicate(timeout=10)
report={'exit':process.returncode,'faults':fail,'stderr':err.decode('utf-8',errors='replace')[-500:]}
(ROOT/'docs/R_NATIVE_FAULT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report),flush=True)
