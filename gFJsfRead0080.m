function [Head,Data]=gFJsfRead0080(JsfHead,ChN,SubSys)
%Read [Head,Data] from JsfHead.fName (*.jsf) file for Message Type 0080 (SBP or SSS Data Message). One sample in packet; take into account Most Significant Bits.
%function [Head,Data]=gFJsfRead0080(JsfHead,ChN,SubSys), where
%JsfHead - Xsf Header structure;
%Head - Header structure;
%ChN - channel number;
%SubSys - subsystem number;
%Data - Data Body for sonar channel number ChN, subsystem number SubSys.
%Head include the next addition fields: HChannelMulti, HSubsystem, HMessageNum.
%The Sonar Data Message consists of a single ping (receiver sounding period) of data for a single channel (such as Port Side Low Frequency Side-Scan).
%Standard sidescan sub-systems have two channels of data, port and starboard.  Standard sub-bottom sub-systems have a single channel of data.
%Data files with higher channel counts exist.  Which fields have data present depends on the system used and data acquisition procedures.
%Example: [Head,Data]=gFJsfRead0080(JsfHead,20,1);

[fId, mes]=fopen(JsfHead.fName,'r');if ~isempty(mes), error(['gFJsfRead0080: ' mes]);end;
LHead=(JsfHead.HMessageType==80)&(JsfHead.HChannelMulti==ChN)&(JsfHead.HSubsystem==SubSys);LenHead=sum(LHead);nHead=find(LHead);
%===Begin Head Allocate for Message Type 80
Head=struct('HChannelMulti',ChN,'HSubsystem',SubSys,'HMessageNum',zeros(1,LenHead),...
    'PingTime',zeros(1,LenHead),'StartingDepth',zeros(1,LenHead),'PingNumber',zeros(1,LenHead),'Reserved1',zeros(2,LenHead),'MSB',zeros(1,LenHead),'Reserved2',zeros(5,LenHead),...
    'IdCode',zeros(1,LenHead),'ValidityFlag',zeros(1,LenHead),'Reserved3',zeros(1,LenHead),'DataFormat',zeros(1,LenHead),'FishAft',zeros(1,LenHead),'FishStb',zeros(1,LenHead),...
    'Reserved4',zeros(2,LenHead),'KilometerPipe',zeros(1,LenHead),'Reserved5',zeros(16,LenHead),'X',zeros(1,LenHead),'Y',zeros(1,LenHead),...
    'CoordinateUnits',zeros(1,LenHead),'AnnotationString',char(zeros(24,LenHead)),'NumberDataSamples',zeros(1,LenHead),'SamplingInterval',zeros(1,LenHead),'GainFactorAdc',zeros(1,LenHead),...
    'PulsePower',zeros(1,LenHead),'Reserved6',zeros(1,LenHead),'ChirpStartingFrequency',zeros(1,LenHead),'ChirpEndingFrequency',zeros(1,LenHead),'SweepLength',zeros(1,LenHead),...
    'Pressure',zeros(1,LenHead),'Depth',zeros(1,LenHead),'SampleFreq',zeros(1,LenHead),'OutgoingPulseId',zeros(1,LenHead),'Altitude',zeros(1,LenHead),'Reserved7',zeros(2,LenHead),...
    'Year',zeros(1,LenHead),'Day',zeros(1,LenHead),'Hour',zeros(1,LenHead),'Minute',zeros(1,LenHead),'Second',zeros(1,LenHead),'TimeBasis',zeros(1,LenHead),...
    'WeightingFactor',zeros(1,LenHead),'NumberOfPulses',zeros(1,LenHead),'CompassHeading',zeros(1,LenHead),'Pitch',zeros(1,LenHead),'Roll',zeros(1,LenHead),'TowElectronicsTemperature',zeros(1,LenHead),...
    'Reserved8',zeros(1,LenHead),'TriggerSource',zeros(1,LenHead),'MarkNumber',zeros(1,LenHead),...
    'NmeaHour',zeros(1,LenHead),'NmeaMinutes',zeros(1,LenHead),'NmeaSeconds',zeros(1,LenHead),'NmeaCourse',zeros(1,LenHead),'NmeaSpeed',zeros(1,LenHead),'NmeaDay',zeros(1,LenHead),'NmeaYear',zeros(1,LenHead),...
    'MillisecondsToday',zeros(1,LenHead),'MaximumAbsoluteValueADC',zeros(1,LenHead),'Reserved9',zeros(1,LenHead),'Reserved10',zeros(1,LenHead),'SoftwareVersionNumber',char(zeros(6,LenHead)),...
    'InitialSphericalCorrectionFactor',zeros(1,LenHead),'PacketNumber',zeros(1,LenHead),'DecimationFactor',zeros(1,LenHead),'DecimationFactorAfterFFT',zeros(1,LenHead),...
    'WaterTemperature',zeros(1,LenHead),'Layback',zeros(1,LenHead),'Reserved11',zeros(1,LenHead),'CableOut',zeros(1,LenHead),'Reserved12',zeros(1,LenHead));
%===Begin Head Allocate for Message Type 80
df=0;
for m=1:LenHead,
    %===Begin Head Read for Message Type 80
    fseek(fId,JsfHead.RSeek(nHead(m))-df,'cof');
    Head.HMessageNum(m)=nHead(m);
    Head.PingTime(m)=fread(fId,1,'int32'); %Ping Time in seconds [since the start of time based on time() function] (1/1/1970) (added in protocol version 8)
    Head.StartingDepth(m)=fread(fId,1,'uint32'); %Starting Depth (window offset) in samples  - usually zero
    Head.PingNumber(m)=fread(fId,1,'uint32'); %Ping Number (increments with each ping)
    Head.Reserved1(:,m)=fread(fId,2,'int16'); %Reserved – Do not use
    Head.MSB(m)=fread(fId,1,'uint16'); %MSBs – Most Significant Bits – High order bits to extend 16 bit unsigned short values to 20 bits.  The 4 MSB bits become the most significant portion of the new 20 bit value. Bits   0 -   3   – start frequency Bits   4 -   7   – end frequency Bits   8 – 11  – samples in this packet Bits 12 – 15  – reserved (added in protocol version 10)(see description below)
    Head.Reserved2(:,m)=fread(fId,5,'int16'); %Reserved – Do not use
    Head.IdCode(m)=fread(fId,1,'int16'); %ID Code (always 1) 1 = Seismic Data
    Head.ValidityFlag(m)=fread(fId,1,'uint16'); %Validity Flag; Validity flags bitmap. Bit0: Lat Lon or XY valid; Bit1: Course valid; Bit2: Speed valid; Bit3: Heading valid; Bit4: Pressure valid; Bit5: Pitch roll valid; Bit6: Altitude valid; Bit7: Reserved; Bit8: Water temperature valid; Bit9: Depth valid; Bit10: Annotation valid; Bit11: Cable counter valid; Bit12: KP valid; Bit13: Position interpolated 
    Head.Reserved3(m)=fread(fId,1,'uint16'); %Reserved – Do not use
    Head.DataFormat(m)=fread(fId,1,'int16'); %Data Format; 0 = 1 short per sample  - Envelope Data; 1 = 2 shorts per sample - Analytic Signal Data, (Real, Imaginary); 2 = 1 short per sample - Raw Data, Prior to Matched Filter; 3 = 1 short per sample - Real portion of Analytic Signal Data; 4 = 1 short per sample - Pixel Data / CEROS Data
    Head.FishAft(m)=fread(fId,1,'int16'); %Distance from Antenna to Tow point in Centimeters, Aft + (Fish Aft = +)
    Head.FishStb(m)=fread(fId,1,'int16'); %Distance from Antenna to Tow Point in Centimeters, Starboard + (Fish to Starboard = +)
    Head.Reserved4(:,m)=fread(fId,2,'int16'); %Reserved – Do not use
    %Navigation data
    Head.KilometerPipe(m)=fread(fId,1,'float32'); %Kilometer of pipe (see bytes 30-31)
    Head.Reserved5(:,m)=fread(fId,16,'int16'); %Reserved – Do not use
    Head.X(m)=fread(fId,1,'int32'); %X in millimeters or decimeters or Longitude in Minutes of Arc / 10000 (see bytes 30-31 and 88-89)
    Head.Y(m)=fread(fId,1,'int32'); %Y in millimeters or decimeters or Latitude in 0.0001 Minutes of Arc (see bytes 30-31 and 88-89)
    Head.CoordinateUnits(m)=fread(fId,1,'int16'); %Coordinate Units: 1 = X, Y in millimeters; 2 = Longitude, Latitude in minutes of arc times 10-4; 3 = X, Y in decimeters
    %Pulse Information
    Head.AnnotationString(:,m)=fread(fId,24,'*char'); %Annotation String (ASCII Data)
    Head.NumberDataSamples(m)=fread(fId,1,'uint16'); %Number of data samples in this packet. See bytes 16 – 17 for MSB information Note: Very large sample sizes require multiple packets
    MSB=uint32(Head.MSB(m));MSB=bitand(MSB,3840);bitshift(MSB,8);Head.NumberDataSamples(m)=Head.NumberDataSamples(m)+MSB;
    Head.SamplingInterval(m)=fread(fId,1,'uint32'); %Sampling Interval in Nanoseconds
    Head.GainFactorAdc(m)=fread(fId,1,'uint16'); %Gain Factor of ADC
    Head.PulsePower(m)=fread(fId,1,'int16'); %User Transmit Level Setting (0 – 100) percen
    Head.Reserved6(m)=fread(fId,1,'int16'); %Reserved – Do not use
    Head.ChirpStartingFrequency(m)=fread(fId,1,'uint16'); %Transmit pulse starting frequency in decahertz (daHz) (units of 10Hz) See bytes 17 – 18 for MSB information
    MSB=uint32(Head.MSB(m));MSB=bitand(MSB,15);bitshift(MSB,16);Head.ChirpStartingFrequency(m)=Head.ChirpStartingFrequency(m)+MSB;
    Head.ChirpEndingFrequency(m)=fread(fId,1,'uint16'); %Transmit pulse ending frequency in decahertz (daHz)(units of 10Hz) See bytes 16 – 17 for MSB information
    MSB=uint32(Head.MSB(m));MSB=bitand(MSB,240);bitshift(MSB,12);Head.ChirpEndingFrequency(m)=Head.ChirpEndingFrequency(m)+MSB;
    Head.SweepLength(m)=fread(fId,1,'uint16'); %Sweep Length in milliseconds
    Head.Pressure(m)=fread(fId,1,'int32'); %Pressure in milliPSI  (1 unit = 1/1000 PSI) (see bytes 30-31)
    Head.Depth(m)=fread(fId,1,'int32'); %Depth in millimeters (if not = 0)  (see bytes 30-31)
    Head.SampleFreq(m)=fread(fId,1,'uint16'); %For all data types EXCEPT RAW (Data Format = 2) this is the Sampling  Frequency of the data. For RAW data, this is one-half the Sample Frequency of the data (Fs/2).  All values are modulo 65536. Use this in conjunction with the Sample interval (Bytes 114-115) to calculate correct sample rate
    Head.OutgoingPulseId(m)=fread(fId,1,'uint16'); %Outgoing pulse identifier
    Head.Altitude(m)=fread(fId,1,'int32'); %Altitude in millimeters (If bottom tracking valid) 0 implies not filled (see bytes 30-31)
    Head.Reserved7(:,m)=fread(fId,2,'int32'); %Reserved – Do not use
    %CPU Time
    Head.Year(m)=fread(fId,1,'int16'); %Year (e.g. 2009) (see Bytes 0-3) (Should not be used)
    Head.Day(m)=fread(fId,1,'int16'); %Day (1 – 366) (Should not be used)
    Head.Hour(m)=fread(fId,1,'int16'); %Hour (see Bytes 200-203) (Should not be used)
    Head.Minute(m)=fread(fId,1,'int16'); %Minute (Should not be used)
    Head.Second(m)=fread(fId,1,'int16'); %Second (should not be used)
    Head.TimeBasis(m)=fread(fId,1,'int16'); %Time Basis (always 3)
    %Weighting Factor
    Head.WeightingFactor(m)=fread(fId,1,'int16'); %Weighting Factor N   (Signed Value!) Defined as 2^-N
    Head.NumberOfPulses(m)=fread(fId,1,'int16'); %Number of pulses in the water 
    %Orientation Sensor Data
    Head.CompassHeading(m)=fread(fId,1,'uint16'); %Compass Heading (0 to 360 ) in units of 1/100 degree (see bytes 30-31)
    Head.Pitch(m)=fread(fId,1,'int16'); %Pitch:  Scale by 180 / 32768 to get degrees, + = bow up (see bytes 30-31)
    Head.Roll(m)=fread(fId,1,'int16'); %Roll:  Scale by 180 / 32768 to get degrees, + = port up (see bytes 30-31)
    Head.TowElectronicsTemperature(m)=fread(fId,1,'int16'); %Tow fish electronics Temperature, in unit of 1/10th degree C
    %Miscellaneous Data
    Head.Reserved8(m)=fread(fId,1,'int16'); %Reserved – Do not use
    Head.TriggerSource(m)=fread(fId,1,'int16'); %Trigger Source: 0 = Internal; 1 = External; 2 = Coupled
    Head.MarkNumber(m)=fread(fId,1,'uint16'); %Mark Number 0=No Mark
    %NMEA Navigation Data
    Head.NmeaHour(m)=fread(fId,1,'int16'); %Hour (0 – 23)
    Head.NmeaMinutes(m)=fread(fId,1,'int16'); %Minutes (0 – 59)
    Head.NmeaSeconds(m)=fread(fId,1,'int16'); %Seconds (0 – 59)
    Head.NmeaCourse(m)=fread(fId,1,'int16'); %Course
    Head.NmeaSpeed(m)=fread(fId,1,'int16'); %Speed
    Head.NmeaDay(m)=fread(fId,1,'int16'); %Day (1 – 366)
    Head.NmeaYear(m)=fread(fId,1,'int16'); %Year
    %Other Miscellaneous Data
    Head.MillisecondsToday(m)=fread(fId,1,'uint32'); %Milliseconds today (since midnight) (use in conjunction with Year / Day to get time of Ping)
    Head.MaximumAbsoluteValueADC(m)=fread(fId,1,'uint16'); %Maximum Absolute Value of ADC samples in this packet
    Head.Reserved9(m)=fread(fId,1,'int16'); %Reserved – Do not use
    Head.Reserved10(m)=fread(fId,1,'int16'); %Reserved – Do not use
    Head.SoftwareVersionNumber(:,m)=fread(fId,6,'*char'); %Sonar Software Version Number - ASCII
    Head.InitialSphericalCorrectionFactor(m)=fread(fId,1,'int32'); % Initial Spherical Correction Factor (Useful for multi-ping / deep application) * 100 
    Head.PacketNumber(m)=fread(fId,1,'uint16'); %Packet Number Each ping starts with packet 1
    Head.DecimationFactor(m)=fread(fId,1,'int16'); %100 times the A/D Decimation Factor. Data is normally sampled at a high Rate.  Digital filters are applied to precisely limit the signal bandwidth.
    Head.DecimationFactorAfterFFT(m)=fread(fId,1,'int16'); %Decimation Factor after the FFT
    Head.WaterTemperature(m)=fread(fId,1,'int16'); %Water Temperature in units of 1/10 degree C (see bytes 30-31)
    Head.Layback(m)=fread(fId,1,'float32'); %Layback in meters
    Head.Reserved11(m)=fread(fId,1,'int32'); %Reserved – Do not use
    Head.CableOut(m)=fread(fId,1,'uint16'); %Cable Out in decimeters (see bytes 30-31) 
    Head.Reserved12(m)=fread(fId,1,'uint16'); %Reserved – Do not use
    %===End Head Read for Message Type 80
    df=ftell(fId);
end;
%===Begin Data Allocate for Message Type 80
if all(Head.DataFormat~=1), Data=zeros(max(Head.NumberDataSamples),LenHead);
else Data=complex(zeros(max(Head.NumberDataSamples),LenHead));
end;
%===End Data Allocate for Message Type 80
fseek(fId,0,'bof');df=0;
for m=1:LenHead,
    if ~mod(m,5000), disp(['Trace: ',num2str(m)]);end;
    %===Begin Data Read for Message Type 80
    fseek(fId,JsfHead.RSeek(nHead(m))+240-df,'cof');
    if Head.DataFormat(m)~=1,
        Data(1:Head.NumberDataSamples(m),m)=fread(fId,Head.NumberDataSamples(m),'int16');
    else
        tmp=fread(fId,Head.NumberDataSamples(m).*2,'int16');
        tmp=reshape(tmp,2,length(tmp)./2);Data(1:Head.NumberDataSamples(m),m)=complex(tmp(1,:),tmp(2,:));
    end;
    %===End Data Read for Message Type 80
    df=ftell(fId);
end;
fclose(fId);

%mail@ge0mlib.ru 01/08/2016