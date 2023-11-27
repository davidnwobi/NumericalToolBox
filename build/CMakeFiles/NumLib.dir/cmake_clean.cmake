file(REMOVE_RECURSE
  "libNumLib.a"
  "libNumLib.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/NumLib.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
