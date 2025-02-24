<script>
	import { config } from '$lib/stores';
	import { get } from 'svelte/store';

	let epsilonStabilityExponent = -7;
	let epsilonKernelExponent = -6;
	let finiteDiffHExponent = -4;
	let constraintToleranceExponent = -7;

	function updateConfig() {
		config.set({
			...get(config),
			epsilonStability: Math.pow(10, epsilonStabilityExponent),
			epsilonKernel: Math.pow(10, epsilonKernelExponent),
			finiteDiffH: Math.pow(10, finiteDiffHExponent),
			constraintTolerance: Math.pow(10, constraintToleranceExponent),
			tauInitial: get(config).tauInitial,
			aConst: get(config).aConst,
			bConst: get(config).bConst,
			maxLineSearch: get(config).maxLineSearch
		});
	}
</script>

<div class="controls" style="display: flex; flex-direction: column; gap: 10px;">
	<label for="epsilonStabilityExponent">
		Epsilon Stability Exponent:
		<input
			type="range"
			id="epsilonStabilityExponent"
			min="-10"
			max="-5"
			step="1"
			bind:value={epsilonStabilityExponent}
			on:input={updateConfig}
		/>
	</label>
	<p>Epsilon Stability: {Math.pow(10, epsilonStabilityExponent).toExponential(2)}</p>

	<label for="epsilonKernelExponent">
		Epsilon Kernel Exponent:
		<input
			type="range"
			id="epsilonKernelExponent"
			min="-10"
			max="-5"
			step="1"
			bind:value={epsilonKernelExponent}
			on:input={updateConfig}
		/>
	</label>
	<p>Epsilon Kernel: {Math.pow(10, epsilonKernelExponent).toExponential(2)}</p>

	<label for="finiteDiffHExponent">
		Finite Diff H Exponent:
		<input
			type="range"
			id="finiteDiffHExponent"
			min="-6"
			max="-3"
			step="1"
			bind:value={finiteDiffHExponent}
			on:input={updateConfig}
		/>
	</label>
	<p>Finite Diff H: {Math.pow(10, finiteDiffHExponent).toExponential(2)}</p>

	<label for="constraintToleranceExponent">
		Constraint Tolerance Exponent:
		<input
			type="range"
			id="constraintToleranceExponent"
			min="-10"
			max="-5"
			step="1"
			bind:value={constraintToleranceExponent}
			on:input={updateConfig}
		/>
	</label>
	<p>Constraint Tolerance: {Math.pow(10, constraintToleranceExponent).toExponential(2)}</p>

	<label for="tauInitial">
		Tau Initial:
		<input
			type="range"
			id="tauInitial"
			min="0.01"
			max="1.0"
			step="0.01"
			bind:value={$config.tauInitial}
			on:input={updateConfig}
		/>
	</label>
	<p>Tau Initial: {$config.tauInitial.toFixed(2)}</p>

	<label for="aConst">
		Armijo Constant:
		<input
			type="range"
			id="aConst"
			min="0.001"
			max="0.5"
			step="0.001"
			bind:value={$config.aConst}
			on:input={updateConfig}
		/>
	</label>
	<p>Armijo Constant: {$config.aConst.toFixed(3)}</p>

	<label for="bConst">
		Backtracking Factor:
		<input
			type="range"
			id="bConst"
			min="0.1"
			max="0.9"
			step="0.01"
			bind:value={$config.bConst}
			on:input={updateConfig}
		/>
	</label>
	<p>Backtracking Factor: {$config.bConst.toFixed(2)}</p>

	<label for="maxLineSearch">
		Max Line Search:
		<input
			type="range"
			id="maxLineSearch"
			min="10"
			max="100"
			step="1"
			bind:value={$config.maxLineSearch}
			on:input={updateConfig}
		/>
	</label>
	<p>Max Line Search: {$config.maxLineSearch}</p>
</div>

<style>
	input[type='range'] {
		width: 200px;
	}
</style>
