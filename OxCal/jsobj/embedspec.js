var s,ss;
spec=new itemSpec("data","List","Object",true);
s=spec.appendChild("group","Group","Object");
s.appendChild("page","Page","Text",true);
s.appendChild("menu","Description","Text",false);
ss=s.appendChild("sub","Sub-group","Array",true);
ss.appendChild("page","Page","Text",true);
ss.appendChild("menu","Description","Text",false);
ss.expand=false;
ss=ss.appendChild("sub","Sub-group","Array",true);
ss.appendChild("page","Page","Text",true);
ss.appendChild("menu","Description","Text",false);
ss.expand=false;
