function signa(waveform,filename,scale)

%
%  signa(waveform, filename [,scale]);
%
%  writes the waveform out as short integers with the low
%  bit masked off.
%
%  Inputs:
%    waveform --  vector, may be complex
%    filename --  string, if wavefrom is complex '.r' and '.i' are appended,
%                    and two files are written.
%    scale    --  optional scale.  If unspecified, the waveform is scaled to
%                   full scale integer 32766.  If specified, the output is
%                   waveform*scale*32766
%  

%
%  Written by John Pauly, Dec. 5, 1994
%  (c) Leland Stanford Jr. University
%

wmax = hex2dec('7ffe');

% if no scale is specified, use as much dynamic range as possible
if nargin == 2,
  scale = 1/max(max(abs(real(waveform)),abs(imag(waveform))));
end;

% scale up to fit in a short integer
waveform = waveform*scale*wmax;

% mask off low bit, since it would be an EOS otherwise
waveform = 2*round(waveform/2);

% if the imaginary component is zero, supress it
if sum(abs(imag(waveform))) == 0,
  waveform = real(waveform);
end;

if isreal(waveform),
  fip = fopen(filename,'wb','b');
  if fip == -1,
    disp(sprintf('Error opening %s for write',filename));
    return;
  end;
  fwrite(fip,waveform,'short');
else
  fip = fopen([filename,'.r'],'wb','b');
  if fip == -1,
    disp(sprintf('Error opening %s for write',[filename,'.r']));
    return;
  end;
  fwrite(fip,real(waveform),'short');
  fclose(fip);
  fip = fopen([filename,'.i'],'wb','b');
  if fip == -1,
    disp(sprintf('Error opening %s for write',[filename,'.i']));
    return;
  end;
  fwrite(fip,imag(waveform),'short');
  fclose(fip);
end;
