echo "Launching Tests for Synaptic Plasticity"
rm *.txt
rm -rf figures_plasticity/*

#### CHECK MODELS #####

python3 test_general.py  &>TestGeneralLog.txt
if [ $? = 0 ]; then
  echo "test_general.py SUCCESS"
else
  echo "test_general.py FAIL"
fi

python3 test_silent_GrC.py  &>TestSilentGrCLog.txt
if [ $? = 0 ]; then
  echo "test_silent_GrC.py SUCCESS"
else
  echo "test_silent_GrC.py FAIL: weight is updated without presynaptic spikes"
fi

python3 test_simple_Gr_spike.py  &>TestSimpleGrCSpikeLog.txt
if [ $? = 0 ]; then
  echo "test_simple_GrC_spike.py SUCCESS"
else
  echo "test_simple_GrC_spike.py FAIL: LTP updates the weight"
fi

python3 test_complex_Gr_spike.py  &>TestComplexGrCSpikeLog.txt
if [ $? = 0 ]; then
  echo "test_complex_Gr_spike.py SUCCESS"
else
  echo "test_complex_Gr_spike.py FAIL: LTD with no presynaptic spike"
fi

python3 test_offset_encoding.py  &>TestOffsetEncodingLog.txt
if [ $? = 0 ]; then
  echo "test_offset_encoding.py SUCCESS"
else
  echo "test_offset_encoding.py FAIL: LTD with no cf spikes"
fi

python3 test_Aplus.py  &>TestAplus.txt
if [ $? = 0 ]; then
  echo "test_Aplus.py SUCCESS"
else
  echo "test_Aplus.py FAIL: Aplus not handled correctly"
fi

python3 test_Aminus.py  &>TestAminus.txt
if [ $? = 0 ]; then
  echo "test_Aminus.py SUCCESS"
else
  echo "test_Aminus.py FAIL: Aminus not handled correctly"
fi

python3 test_null_parameters.py  &>TestNullParameters.txt
if [ $? = 0 ]; then
  echo "test_null_parameters.py SUCCESS"
else
  echo "test_null_parameters.py FAIL: plasticity with null Aplus and Aminus"
fi

python3 test_sinexp_kernel.py  &>TestSinExpKernel.txt
if [ $? = 0 ]; then
  echo "test_sinexp_kernel.py SUCCESS"
else
  echo "test_sinexp_kernel.py FAIL"
fi
