
function(get_cpu_vendor CPU_DESCRIPTION CPU_VENDOR)
  # Function's primary purpose is to identify the CPU Vendor
  # on the Host system, so the correct backend for BLAS/LAPACK
  # and/or FFT fwill be used. Possible values are:
  # AQMD, Intel or Unknown.
  message(STATUS "Attempting to figure out your CPU...")
  #unset(CPU_VENDOR)
  set(_cpus_ "Intel" "AMD")
  foreach(_regex_ IN LISTS _cpus_)
    string(REGEX MATCH ${_regex_} _cpu_vendor_ "${CPU_DESCRIPTION}")
    if(_cpu_vendor_ MATCHES "Intel")
      set(_vendor_name_ "Intel")
    elseif(_cpu_vendor_ MATCHES "AMD")
      set(_vendor_name_ "AMD")
    endif()
  endforeach()

  if(NOT _vendor_name_)
    set(_vendor_name_ "Unknown")
  endif()
  
  message(STATUS "Your CPU is ${_vendor_name_}")
  # Returning the Vendor's name
  set(${CPU_VENDOR} "${_vendor_name_}" PARENT_SCOPE)
endfunction()
