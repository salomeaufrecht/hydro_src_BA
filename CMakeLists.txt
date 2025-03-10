#add_executable("firstexample" firstexample.cc)
#target_link_dune_default_libraries(firstexample)
#target_link_libraries(firstexample PUBLIC ${GDAL_LIBRARIES})
#target_include_directories(firstexample PUBLIC ${GDAL_INCLUDE_DIRS})
#
#add_executable("secondexample" secondexample.cc)
#target_link_dune_default_libraries(secondexample)
#target_link_libraries(secondexample PUBLIC ${GDAL_LIBRARIES})
#target_include_directories(secondexample PUBLIC ${GDAL_INCLUDE_DIRS})
#
#add_executable("thirdexample" thirdexample.cc)
#target_link_dune_default_libraries(thirdexample)
#target_link_libraries(thirdexample PUBLIC ${GDAL_LIBRARIES})
#target_include_directories(thirdexample PUBLIC ${GDAL_INCLUDE_DIRS})
#
#add_executable("fourthexample" fourthexample.cc)
#target_link_dune_default_libraries(fourthexample)
#target_link_libraries(fourthexample PUBLIC GDAL::GDAL)
#
#add_executable("readtiff" readtiff.cc)
#target_link_dune_default_libraries(readtiff)
#target_link_libraries(readtiff PUBLIC GDAL::GDAL)
#target_link_libraries(readtiff PUBLIC NetCDF::NetCDF)


#Added
#add_executable("test_grid_functions" test_grid_functions.cc)
#target_link_dune_default_libraries(test_grid_functions)
#target_link_libraries(test_grid_functions PUBLIC GDAL::GDAL)
#target_link_libraries(test_grid_functions PUBLIC NetCDF::NetCDF)

#add_executable("test02" test02.cc)
#target_link_dune_default_libraries(test02)
#target_link_libraries(test02 PUBLIC GDAL::GDAL)
#target_link_libraries(test02 PUBLIC NetCDF::NetCDF)

#add_executable("test00" test00.cc)
#target_link_dune_default_libraries(test00)
#target_link_libraries(test00 PUBLIC GDAL::GDAL)
#target_link_libraries(test00 PUBLIC NetCDF::NetCDF)

#add_executable("test03" test03.cc)
#target_link_dune_default_libraries(test03)
#target_link_libraries(test03 PUBLIC GDAL::GDAL)
#target_link_libraries(test03 PUBLIC NetCDF::NetCDF)

#add_executable("test04" test04.cc)
#target_link_dune_default_libraries(test04)
#target_link_libraries(test04 PUBLIC GDAL::GDAL)
#target_link_libraries(test04 PUBLIC NetCDF::NetCDF)

#add_executable("test05_borderInTriangle" test05_borderInTriangle.cc)
#target_link_dune_default_libraries(test05_borderInTriangle)
#target_link_libraries(test05_borderInTriangle PUBLIC GDAL::GDAL)
#target_link_libraries(test05_borderInTriangle PUBLIC NetCDF::NetCDF)

#add_executable("test06_borderInTriangle" test06_borderInTriangle.cc)
#target_link_dune_default_libraries(test06_borderInTriangle)
#target_link_libraries(test06_borderInTriangle PUBLIC GDAL::GDAL)
#target_link_libraries(test06_borderInTriangle PUBLIC NetCDF::NetCDF)

#add_executable("test07_diagonal" test07_diagonal.cc)
#target_link_dune_default_libraries(test07_diagonal)
#target_link_libraries(test07_diagonal PUBLIC GDAL::GDAL)
#target_link_libraries(test07_diagonal PUBLIC NetCDF::NetCDF)

#add_executable("test08_overallFlow" test08_overallFlow.cc)
#target_link_dune_default_libraries(test08_overallFlow)
#target_link_libraries(test08_overallFlow PUBLIC GDAL::GDAL)
#target_link_libraries(test08_overallFlow PUBLIC NetCDF::NetCDF)

#add_executable("test09_headerData" test09_headerData.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(test09_headerData)
#target_link_libraries(test09_headerData PUBLIC GDAL::GDAL)
#target_link_libraries(test09_headerData PUBLIC NetCDF::NetCDF)

#add_executable("test10_quantityTriangles" test10_quantityTriangles.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(test10_quantityTriangles)
#target_link_libraries(test10_quantityTriangles PUBLIC GDAL::GDAL)
#target_link_libraries(test10_quantityTriangles PUBLIC NetCDF::NetCDF)
#
#add_executable("test11_seperatetesting" test11_seperatetesting.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(test11_seperatetesting)
#target_link_libraries(test11_seperatetesting PUBLIC GDAL::GDAL)
#target_link_libraries(test11_seperatetesting PUBLIC NetCDF::NetCDF)


#add_executable("test06" test06.cc)
#target_link_dune_default_libraries(test06)
#target_link_libraries(test06 PUBLIC GDAL::GDAL)
#target_link_libraries(test06 PUBLIC NetCDF::NetCDF)

#add_executable("detect_rivers01" detect_rivers01.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(detect_rivers01)
#target_link_libraries(detect_rivers01 PUBLIC GDAL::GDAL)
#target_link_libraries(detect_rivers01 PUBLIC NetCDF::NetCDF)

#add_executable("detect_rivers02" detect_rivers02.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(detect_rivers02)
#target_link_libraries(detect_rivers02 PUBLIC GDAL::GDAL)
#target_link_libraries(detect_rivers02 PUBLIC NetCDF::NetCDF)

#add_executable("detectRivers_countFragmetns" detectRivers_countFragmetns.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(detectRivers_countFragmetns)
#target_link_libraries(detectRivers_countFragmetns PUBLIC GDAL::GDAL)
#target_link_libraries(detectRivers_countFragmetns PUBLIC NetCDF::NetCDF)

#add_executable("detect_upwards_rivers" detect_upwards_rivers.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(detect_upwards_rivers)
#target_link_libraries(detect_upwards_rivers PUBLIC GDAL::GDAL)
#target_link_libraries(detect_upwards_rivers PUBLIC NetCDF::NetCDF)

#add_executable("detect_rivers_example" detect_rivers_example.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(detect_rivers_example)
#target_link_libraries(detect_rivers_example PUBLIC GDAL::GDAL)
#target_link_libraries(detect_rivers_example PUBLIC NetCDF::NetCDF)

#add_executable("detectRivers_oneFunction" detectRivers_oneFunction.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(detectRivers_oneFunction)
#target_link_libraries(detectRivers_oneFunction PUBLIC GDAL::GDAL)
#target_link_libraries(detectRivers_oneFunction PUBLIC NetCDF::NetCDF)

#add_executable("detectRivers_example_ausarbeitung" detectRivers_example_ausarbeitung.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(detectRivers_example_ausarbeitung)
#target_link_libraries(detectRivers_example_ausarbeitung PUBLIC GDAL::GDAL)
#target_link_libraries(detectRivers_example_ausarbeitung PUBLIC NetCDF::NetCDF)

add_executable("bachelorarbeit_withOriginalMap" bachelorarbeit_withOriginalMap.cc flowFunctions.cc flowFunctions.hh)
target_link_dune_default_libraries(bachelorarbeit_withOriginalMap)
target_link_libraries(bachelorarbeit_withOriginalMap PUBLIC GDAL::GDAL)
target_link_libraries(bachelorarbeit_withOriginalMap PUBLIC NetCDF::NetCDF)

#add_executable("bachelorarbeit_beispiel" bachelorarbeit_beispiel.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(bachelorarbeit_beispiel)
#target_link_libraries(bachelorarbeit_beispiel PUBLIC GDAL::GDAL)
#target_link_libraries(bachelorarbeit_beispiel PUBLIC NetCDF::NetCDF)

#add_executable("explain_refinement" explain_refinement.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(explain_refinement)
#target_link_libraries(explain_refinement PUBLIC GDAL::GDAL)
#target_link_libraries(explain_refinement PUBLIC NetCDF::NetCDF)

#add_executable("analyse_1x1" analyse_1x1.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(analyse_1x1)
#target_link_libraries(analyse_1x1 PUBLIC GDAL::GDAL)
#target_link_libraries(analyse_1x1 PUBLIC NetCDF::NetCDF)

#add_executable("analyse_1x1_diag" analyse_1x1_diag.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(analyse_1x1_diag)
#target_link_libraries(analyse_1x1_diag PUBLIC GDAL::GDAL)
#target_link_libraries(analyse_1x1_diag PUBLIC NetCDF::NetCDF)

#add_executable("analyse_diag" analyse_diag.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(analyse_diag)
#target_link_libraries(analyse_diag PUBLIC GDAL::GDAL)
#target_link_libraries(analyse_diag PUBLIC NetCDF::NetCDF)

#add_executable("analyse_1x1_fixedWidth" analyse_1x1_fixedWidth.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(analyse_1x1_fixedWidth)
#target_link_libraries(analyse_1x1_fixedWidth PUBLIC GDAL::GDAL)
#target_link_libraries(analyse_1x1_fixedWidth PUBLIC NetCDF::NetCDF)

#add_executable("analyse_3x3" analyse_3x3.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(analyse_3x3)
#target_link_libraries(analyse_3x3 PUBLIC GDAL::GDAL)
#target_link_libraries(analyse_3x3 PUBLIC NetCDF::NetCDF)

#add_executable("runtime_detectFragments" runtime_detectFragments.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(runtime_detectFragments)
#target_link_libraries(runtime_detectFragments PUBLIC GDAL::GDAL)
#target_link_libraries(runtime_detectFragments PUBLIC NetCDF::NetCDF)

#add_executable("runtime_all" runtime_all.cc flowFunctions.cc flowFunctions.hh)
#target_link_dune_default_libraries(runtime_all)
#target_link_libraries(runtime_all PUBLIC GDAL::GDAL)
#target_link_libraries(runtime_all PUBLIC NetCDF::NetCDF)

#add_executable("sphere" sphere.cc)
#target_link_dune_default_libraries(sphere)
#target_link_libraries(sphere PUBLIC GDAL::GDAL)
#target_link_libraries(sphere PUBLIC NetCDF::NetCDF)

#End added



#add_executable("nakhon" nakhon.cc)
#target_link_dune_default_libraries(nakhon)
#target_link_libraries(nakhon PUBLIC GDAL::GDAL)
#target_link_libraries(nakhon PUBLIC NetCDF::NetCDF)
#
#add_executable("bw90mpre" bw90mpre.cc)
#target_link_dune_default_libraries(bw90mpre)
#target_link_libraries(bw90mpre PUBLIC GDAL::GDAL)
#target_link_libraries(bw90mpre PUBLIC NetCDF::NetCDF)
#
#add_executable("stok" stok.cc)
#target_link_dune_default_libraries(stok)
#target_link_libraries(stok PUBLIC GDAL::GDAL)
#target_link_libraries(stok PUBLIC NetCDF::NetCDF)
#
#add_executable("stok_long_term" stok_long_term.cc)
#target_link_dune_default_libraries(stok_long_term)
#target_link_libraries(stok_long_term PUBLIC GDAL::GDAL)
#target_link_libraries(stok_long_term PUBLIC NetCDF::NetCDF)
#
#add_executable("coupledexample" coupledexample.cc)
#target_link_dune_default_libraries(coupledexample)
#target_link_libraries(coupledexample PUBLIC GDAL::GDAL)
#target_link_libraries(coupledexample PUBLIC NetCDF::NetCDF)
#
#add_executable("VCFV_example" VCFV_example.cc)
#target_link_dune_default_libraries(VCFV_example)
#target_link_libraries(VCFV_example PUBLIC ${GDAL_LIBRARIES})
#target_include_directories(VCFV_example PUBLIC ${GDAL_INCLUDE_DIRS})
#
dune_symlink_to_source_files(FILES fourthexample.ini thirdexample.ini secondexample.ini firstexample.ini readtiff.ini nakhon.ini bw90mpre.ini stok.ini coupledexample.ini stok_long_term.ini VCFV_example.ini)
#
#add_executable("dune-hydro" dune-hydro.cc)
#target_link_dune_default_libraries(dune-hydro)
#target_link_libraries(dune-hydro PUBLIC ${GDAL_LIBRARIES})
#target_include_directories(dune-hydro PUBLIC ${GDAL_INCLUDE_DIRS})
#