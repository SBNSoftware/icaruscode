#!/usr/bin/perl -w

sub read_pmt_pos {

  $pmt_x = "0";

  $PMT_pos_file = $_[0];
  open(PMTPOS, $PMT_pos_file) or die("Could not open file $PMT_pos_file.");

  @pmt_pos = ();

  foreach $line (<PMTPOS>) {

    @coord = split(/\s/, $line);

    $string = "' x=\" $pmt_x\" y=\"$coord[0]\" z=\"$coord[1]\"'"; 
    push(@pmt_pos, $string);

  }

  close(PMTPOS);

  return @pmt_pos;
}

@pmt_pos = read_pmt_pos("disposizione_GL.txt");
print join(",\n", @pmt_pos);
