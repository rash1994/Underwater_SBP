function JsfHead=gFJsfReadHeader(fName,flStat)
%Read JsfHead structure from *.jsf file.
%function JsfHead=gFJsfReadHeader(fName,flStat), where
%JsfHead - JsfHead structure;
%fName - the target file name;
%flStat - flag for statistics display (1 or 0).
%JsfHead include the next addition fields: JsfHead.Descript, JsfHead.fName, JsfHead.RSeek, JsfHead.ROnFlag.
%Example: JsfHead=gFJsfReadHeader('c:\temp\1.jsf',1);

Descript.MessageType.Code=[80,82,86,182,426,428,2000,2002,2020,2040,2060,2080,2090,2091,2100,2101,2111,9001,9002,9003]';
Descript.MessageType.Text={'80=Sonar Data Message','82=Side Scan Data Message','86=4400-SAS Processed Data',...
    '182=System Information Message',...
    '426=File Timestamp Message','428=File Padding Message',...
    '2000=Equip Serial Ports Raw Data','2002=NMEA String','2020=Pitch Roll Data','2040=Miscellaneous Analog Sensors','2060=Pressure Sensor Reading','2080=Doppler Velocity Log Data (DVL)','2090=Situation Message','2091=Situation Comprehensive Message (version 2)',...
    '2100=Cable Counter Data Message','2101=Kilometer of Pipe Data','2111=Container Timestamp Message',...
    '9001=Discover-2 General Prefix Message','9002=Discover-2 Situation Data','9003=Discover-2 Acoustic Prefix Message'}';
Descript.Subsystem.Code=[0,20,21,100,101,102]';
Descript.Subsystem.Text={'0=Sub-bottom','20=Single frequency sidescan data or Lower Frequency Dual Side Scan (75 or 120kHz typical)','21=Higher Frequency Dual Side Scan (410kHz typical)','100=Raw serial data','101=Parsed serial data','102=Miscellaneous Analog Sensors (Mess2040)'}';
Descript.ChannelMulti.Code=[0,1]';
Descript.ChannelMulti.Text={'0=Port','1=Starboard'}';
Descript.SystemTypeNumber.Code=[1,2,4,5,6,7,11,14,16,17,18,19,20,21,23,24,25,27,128]';
Descript.SystemTypeNumber.Text={'2xxx Series, Combined Sub-Bottom / Side Scan with SIB Electronics','2xxx Series, Combined Sub-Bottom / Side Scan with FSIC Electronics','4300-MPX(Multi-Ping','3200-XS,Sub-Bottom Profiler wit h AIC Electronics',...
        '4400-SAS, 12-Channel Side Scan','3200-XS, Sub Bottom Profiler with SIB Electronics','4200 Limited Multipulse Dual Frequency Side Scan','3100-P, Sub Bottom Profiler','2xxx Series, Dual Side Scan with SIB Electronics',...
        '4200 Multipulse Dual Frequency Side Scan','4700 Dynamic Focus','4200 Dual Frequency Side Scan','4200 Dual Frequency non Simultaneous Side Scan','2200-MP Combined Sub-Bottom / Dual Frequency Multipulse Side Scan',...
        '4600 Multipulse Bathymetric System','4200 Single Frequency Dynamically Focused Side Scan','4125 Dual Frequency Side Scan','4600 Monopulse Bathymetric Syste','4100, 272 /560A Side Scan'};

[fId, mes]=fopen(fName,'r');if ~isempty(mes), error(['gFJsfReadHeader: ' mes]);end;
%===Begin Calc Num of Records
nRec=0;finfo=dir(fName);fSize=finfo.bytes;
while fSize>ftell(fId),
    nRec=nRec+1;
    face=fread(fId,1,'uint16');if face~=5633, error('Error gFJsfReadHeader: Marker for the Start of JstHeader~=0x1601');end;
    fseek(fId,10,'cof');
    SizeFollowingMessage=fread(fId,1,'uint32'); %Size of following Message in Bytes
    fseek(fId,SizeFollowingMessage,'cof'); %goto block end
end;
%===End Calc Num of Records
%===Begin JsfHeader Record Allocate
JsfHead=struct('Descript',Descript,'fName',fName,'HMarkerForStart',nan(1,nRec),'HVersionOfProtocol',nan(1,nRec),'HSessionIdentifier',nan(1,nRec),'HMessageType',nan(1,nRec),'HCommandType',nan(1,nRec),...
    'HSubsystem',nan(1,nRec),'HChannelMulti',nan(1,nRec),'HSequenceNumber',nan(1,nRec),'HReserved',nan(1,nRec),'HSizeFollowingMessage',nan(1,nRec),...
    'RSeek',nan(1,nRec),'ROnFlag',nan(1,nRec));
%===End JsfHeader Record Allocate
fseek(fId,0,'bof');
for m=1:nRec,
    %===Begin JsfHeader Record Read
    JsfHead.HMarkerForStart(m)=fread(fId,1,'uint16'); %Marker for the Start of Header = 0x1601
    JsfHead.HVersionOfProtocol(m)=fread(fId,1,'uint8'); %Version of Protocol used
    JsfHead.HSessionIdentifier(m)=fread(fId,1,'uint8'); %Session Identifier
    JsfHead.HMessageType(m)=fread(fId,1,'uint16'); %Message Type
    JsfHead.HCommandType(m)=fread(fId,1,'uint8'); %Command Type
    JsfHead.HSubsystem(m)=fread(fId,1,'uint8'); %Subsystem for a Multi-System Device 0=Sub-bottom; 20=75 or 120 kHz Side Scan; 21=410 kHz Side Scan
    JsfHead.HChannelMulti(m)=fread(fId,1,'uint8'); %Channel for a Multi-Channel Subsystem For Side Scan Subsystems; 0 = Port; 1 = Starboard
    JsfHead.HSequenceNumber(m)=fread(fId,1,'uint8'); %Sequence Number
    JsfHead.HReserved(m)=fread(fId,1,'uint16'); %Reserved
    JsfHead.HSizeFollowingMessage(m)=fread(fId,1,'uint32'); %Size of following Message in Bytes
    JsfHead.RSeek(m)=ftell(fId); %Message Seeker
    JsfHead.ROnFlag(m)=1; %On/Off flag
    %===End JsfHeader Record Read
    fseek(fId,JsfHead.HSizeFollowingMessage(m),'cof'); %goto block end
end;
fclose(fId);

%===Begin Statistics display
if flStat,
    a1=sparse(JsfHead.HMessageType+1,ones(1,nRec),ones(1,nRec));a1Mess=find(a1)-1;a1Num=nonzeros(a1);
    for n=1:size(a1Mess,1),
        if any(JsfHead.Descript.MessageType.Code==a1Mess(n)), s=JsfHead.Descript.MessageType.Text{JsfHead.Descript.MessageType.Code==a1Mess(n)};
        else s=[num2str(a1Mess(n)) '=Not Defined'];end;
        fprintf('Mess: %s, Num: %d ',s,a1Num(n));
        L=find(JsfHead.HMessageType==a1Mess(n));L2=size(L,2);
        a2=sparse(JsfHead.HSubsystem(L)+1,ones(1,L2),ones(1,L2));a2Mess=find(a2)-1;a2Num=nonzeros(a2);
        for nn=1:size(a2Mess,1),
            fprintf('[ Subs: %d, Num: %d;',a2Mess(nn),a2Num(nn));
            L=find((JsfHead.HMessageType==a1Mess(n))&(JsfHead.HSubsystem==a2Mess(nn)));L2=size(L,2);
            a3=sparse(JsfHead.HChannelMulti(L)+1,ones(1,L2),ones(1,L2));a3Mess=find(a3)-1;a3Num=nonzeros(a3);
            fprintf(' Chan:');fprintf(' %d',a3Mess);fprintf(', Num:');fprintf(' %d',a3Num);fprintf(' ]');
        end;
        fprintf('\n');
    end;
end;
%===End Statistics display

%mail@ge0mlib.ru 01/08/2016