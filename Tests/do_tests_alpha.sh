echo "Launching Tests for Synaptic Plasticity"
rm *.txt
rm -rf figures_plasticity_alpha_normal_synapse/*

#### CHECK MODELS #####

python3 test_general_alpha.py  &>TestGeneralLog_alpha.txt
if [ $? = 0 ]; then
  echo "test_general_alpha.py SUCCESS"
else
  echo "test_general_alpha.py FAIL"
fi

python3 test_silent_GrC_alpha.py  &>TestSilentGrCLog_alpha.txt
if [ $? = 0 ]; then
  echo "test_silent_GrC_alpha.py SUCCESS"
else
  echo "test_silent_GrC_alpha.py FAIL: weight is updated without presynaptic spikes"
fi

python3 test_simple_Gr_spike_alpha.py  &>TestSimpleGrCSpikeLog_alpha.txt
if [ $? = 0 ]; then
  echo "test_simple_GrC_spike_alpha.py SUCCESS"
else
  echo "test_simple_GrC_spike_alpha.py FAIL: LTD updates the weight"
fi

python3 test_complex_Gr_spike_alpha.py  &>TestComplexGrCSpikeLog_alpha.txt
if [ $? = 0 ]; then
  echo "test_complex_Gr_spike_alpha.py SUCCESS"
else
  echo "test_complex_Gr_spike_alpha.py FAIL: LTP with no presynaptic spike"
fi

python3 test_offset_encoding_alpha.py  &>TestOffsetEncodingLog_alpha.txt
if [ $? = 0 ]; then
  echo "test_offset_encoding_alpha.py SUCCESS"
else
  echo "test_offset_encoding_alpha.py FAIL: LTP with no cf spikes"
fi

python3 test_Aplus_alpha.py  &>TestAplus_alpha.txt
if [ $? = 0 ]; then
  echo "test_Aplus_alpha.py SUCCESS"
else
  echo "test_Aplus_alpha.py FAIL: Aplus not handled correctly"
fi

python3 test_Aminus_alpha.py  &>TestAminus_alpha.txt
if [ $? = 0 ]; then
  echo "test_Aminus_alpha.py SUCCESS"
else
  echo "test_Aminus_alpha.py FAIL: Aminus not handled correctly"
fi

python3 test_null_parameters_alpha.py  &>TestNullParameters_alpha.txt
if [ $? = 0 ]; then
  echo "test_null_parameters_alpha.py SUCCESS"
else
  echo "test_null_parameters_alpha.py FAIL: plasticity with null Aplus and Aminus"
fi

python3 test_alpha_kernel.py  &>TestAlphaKernel.txt
if [ $? = 0 ]; then
  echo "test_alpha_kernel_alpha.py SUCCESS"
else
  echo "test_alpha_kernel.py FAIL"
fi