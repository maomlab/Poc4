#! /bin/awk
# From https://bioinformatics.stackexchange.com/a/21368
# Usage:
#   awk -f m8_to_fasta.awk my_m8_file.m8 > my_fasta_file.fasta

{
  name = "> "
  # start on the second field, foldseek will put the job name in $1
  field = 2
  # Check each field until we see numbers
  while ( $field + 0 != $field ) {
     name = name " " $field
     field = field+1
  }
  # name now contains the name string which is the first line of a fasta entry

  field = NF
  # Now find the matched sequence
  # Sometimes an entry ENDS in a number, so we strip that
  if ( $field + 0 == $field ) { field = field -1 }
  while ( $field + 0 != $field ) {
      field = field - 1
  }
  # field is now pointing at some undocumented number. Decrement once more to get the sequence
  field = field - 1
  # field is now pointing at the sequence
  seq = $field

  #Now I/we want to sanitize the sequence to remove any gaps
  gsub(/-/, "", seq)

  print name
  print seq
}
