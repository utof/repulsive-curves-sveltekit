<script>
	import { createEventDispatcher } from 'svelte';
	import {
		config,
		optimizationConfig,
		initialTotalLength,
		currentUser,
		currentDateTimeUTC
	} from './stores';
	import { GradientMethods } from './optimization';
	import { get } from 'svelte/store';

	const dispatch = createEventDispatcher();

	// Local copies of store values for binding
	// Alpha beta parameters
	let alpha = $optimizationConfig.alpha;
	let beta = $optimizationConfig.beta;
	let linkAlphaBeta = $optimizationConfig.linkAlphaBeta;

	// Optimization settings
	let gradientMethod = $optimizationConfig.gradientMethod;

	// Constraint settings
	let useBarycenterConstraint = $optimizationConfig.constraints.barycenter.enabled;
	let barycenterTarget = [...$optimizationConfig.constraints.barycenter.target];
	let useLengthConstraint = $optimizationConfig.constraints.length.enabled;
	let lengthTarget = $optimizationConfig.constraints.length.absoluteValue;
	let lengthPercentage = $optimizationConfig.constraints.length.percentage;
	let useLengthPercentage = $optimizationConfig.constraints.length.usePercentage;
	let useFullConstraintProjection = $optimizationConfig.useFullConstraintProjection;

	// Step size controls
	let useLineSearch = $optimizationConfig.useLineSearch;
	let precondStepSize = $optimizationConfig.precondStepSize;
	let l2StepSize = $optimizationConfig.l2StepSize;

	// Constraint stabilization
	let barycenterStabilization = $optimizationConfig.barycenterStabilization;
	let lengthStabilization = $optimizationConfig.lengthStabilization;

	// Panel toggle states
	let openPanels = {
		optimize: true,
		constraints: true,
		numerical: false,
		advanced: false
	};

	// Exponent sliders
	let epsilonStabilityExponent = Math.log10($config.epsilonStability);
	let epsilonKernelExponent = Math.log10($config.epsilonKernel);
	let finiteDiffHExponent = Math.log10($config.finiteDiffH);
	let constraintToleranceExponent = Math.log10($config.constraintTolerance);

	// Stores reference to the optimizer
	let optimizer;

	// Export a function to set the optimizer reference
	export function setOptimizer(optimizerInstance) {
		optimizer = optimizerInstance;
		updateOptimizerConfig();
	}

	// Update numerical configuration parameters
	function updateConfig() {
		config.update((current) => ({
			...current,
			epsilonStability: Math.pow(10, epsilonStabilityExponent),
			epsilonKernel: Math.pow(10, epsilonKernelExponent),
			finiteDiffH: Math.pow(10, finiteDiffHExponent),
			constraintTolerance: Math.pow(10, constraintToleranceExponent)
		}));

		dispatch('update');
	}

	// Update alpha/beta settings
	function updateAlphaBeta() {
		if (linkAlphaBeta) {
			beta = 2 * alpha;
		}

		optimizationConfig.update((current) => ({
			...current,
			alpha,
			beta,
			linkAlphaBeta
		}));

		dispatch('alphaBetaChange', { alpha, beta });
		dispatch('update');
	}

	// Update optimization configuration
	function updateOptimizationSettings() {
		optimizationConfig.update((current) => ({
			...current,
			gradientMethod,
			useLineSearch,
			precondStepSize,
			l2StepSize,
			barycenterStabilization,
			lengthStabilization,
			useFullConstraintProjection,
			constraints: {
				barycenter: {
					enabled: useBarycenterConstraint,
					target: barycenterTarget
				},
				length: {
					enabled: useLengthConstraint,
					usePercentage: useLengthPercentage,
					percentage: lengthPercentage,
					absoluteValue: lengthTarget
				}
			}
		}));

		dispatch('update');
	}

	// Update optimizer configuration when settings change
	function updateOptimizerConfig() {
		if (!optimizer) return;

		// Update gradient method
		optimizer.setGradientMethod(gradientMethod);

		// Update constraints
		optimizer.setConstraints({
			barycenter: useBarycenterConstraint,
			barycenterTarget,
			length: useLengthConstraint,
			lengthTarget: useLengthPercentage ? undefined : lengthTarget,
			lengthPercentage: useLengthPercentage ? lengthPercentage : undefined
		});

		// Update constraint projection method
		optimizer.setUseFullConstraintProjection(useFullConstraintProjection);

		// Update any other optimizer settings
		optimizer.updateSettings({
			useLineSearch,
			precondStepSize,
			l2StepSize
		});

		// Notify parent component
		dispatch('update');
	}

	// Toggle panel visibility
	function togglePanel(panel) {
		openPanels[panel] = !openPanels[panel];
	}

	// When optimization config store changes, update local variables
	function updateFromStore() {
		const opt = get(optimizationConfig);
		alpha = opt.alpha;
		beta = opt.beta;
		linkAlphaBeta = opt.linkAlphaBeta;
		gradientMethod = opt.gradientMethod;
		useBarycenterConstraint = opt.constraints.barycenter.enabled;
		barycenterTarget = [...opt.constraints.barycenter.target];
		useLengthConstraint = opt.constraints.length.enabled;
		lengthTarget = opt.constraints.length.absoluteValue;
		lengthPercentage = opt.constraints.length.percentage;
		useLengthPercentage = opt.constraints.length.usePercentage;
		useFullConstraintProjection = opt.useFullConstraintProjection;
		useLineSearch = opt.useLineSearch;
		precondStepSize = opt.precondStepSize;
		l2StepSize = opt.l2StepSize;
		barycenterStabilization = opt.barycenterStabilization;
		lengthStabilization = opt.lengthStabilization;
	}

	// Subscribe to store changes
	const unsubscribe = optimizationConfig.subscribe(() => {
		updateFromStore();
	});

	// Clean up subscription
	import { onDestroy } from 'svelte';
	onDestroy(unsubscribe);

	// Apply changes when optimization settings are updated
	$: {
		// This reactive block triggers when any optimization setting changes
		if ($optimizationConfig) {
			// This only updates the UI without triggering loops
			if (optimizer) {
				updateOptimizerConfig();
			}
		}
	}

	// Format a number with commas
	function formatNumber(num) {
		return num.toLocaleString(undefined, { maximumFractionDigits: 2 });
	}

	// Convert exponents to readable scientific notation
	function formatExponent(exponent) {
		return Math.pow(10, exponent).toExponential(2);
	}

	// Format the current date for display
	function formatDate(dateString) {
		return new Date(dateString).toLocaleString();
	}
</script>

<div class="controls-container">
	<div class="controls-header">
		<h3>Optimization Controls</h3>
		<div class="user-info">
			<span>User: {$currentUser}</span>
			<span title={formatDate($currentDateTimeUTC)}>
				{$currentDateTimeUTC.split(' ')[0]}
			</span>
		</div>
	</div>

	<!-- OPTIMIZATION PARAMETERS SECTION -->
	<div class="control-panel">
		<div class="panel-header" on:click={() => togglePanel('optimize')}>
			<h4>Optimization Parameters</h4>
			<span class="toggle-icon">{openPanels.optimize ? '▼' : '►'}</span>
		</div>

		{#if openPanels.optimize}
			<div class="panel-content">
				<!-- Alpha & Beta Controls -->
				<div class="control-group">
					<label>
						Alpha:
						<input
							type="number"
							bind:value={alpha}
							on:change={updateAlphaBeta}
							min="1"
							max="10"
							step="0.5"
						/>
					</label>

					<label>
						Beta:
						<input
							type="number"
							bind:value={beta}
							on:change={updateAlphaBeta}
							min="2"
							max="20"
							step="0.5"
							disabled={linkAlphaBeta}
						/>
					</label>

					<label class="checkbox-label">
						<input type="checkbox" bind:checked={linkAlphaBeta} on:change={updateAlphaBeta} />
						Set β = 2α (paper recommendation)
					</label>
				</div>

				<!-- Gradient Method -->
				<div class="control-group">
					<h5>Gradient Method</h5>
					<label class="radio-label">
						<input
							type="radio"
							bind:group={gradientMethod}
							value={GradientMethods.PRECONDITIONED}
							on:change={updateOptimizationSettings}
						/>
						Fractional Sobolev (Paper Method)
					</label>
					<label class="radio-label">
						<input
							type="radio"
							bind:group={gradientMethod}
							value={GradientMethods.L2}
							on:change={updateOptimizationSettings}
						/>
						Standard L2 Gradient
					</label>
				</div>

				<!-- Step Size Controls -->
				<div class="control-group">
					<h5>Step Size Control</h5>
					<label class="checkbox-label">
						<input
							type="checkbox"
							bind:checked={useLineSearch}
							on:change={updateOptimizationSettings}
						/>
						Use Line Search
					</label>

					{#if !useLineSearch}
						<div class="slider-group">
							<label>
								Preconditioned Step Size:
								<div class="slider-with-value">
									<input
										type="range"
										bind:value={precondStepSize}
										on:change={updateOptimizationSettings}
										min="0.000001"
										max="0.0001"
										step="0.000001"
									/>
									<span>{precondStepSize.toExponential(2)}</span>
								</div>
							</label>

							<label>
								L2 Step Size:
								<div class="slider-with-value">
									<input
										type="range"
										bind:value={l2StepSize}
										on:change={updateOptimizationSettings}
										min="10000"
										max="1000000"
										step="10000"
									/>
									<span>{l2StepSize.toExponential(2)}</span>
								</div>
							</label>
						</div>
					{/if}
				</div>
			</div>
		{/if}
	</div>

	<!-- CONSTRAINTS SECTION -->
	<div class="control-panel">
		<div class="panel-header" on:click={() => togglePanel('constraints')}>
			<h4>Constraints</h4>
			<span class="toggle-icon">{openPanels.constraints ? '▼' : '►'}</span>
		</div>

		{#if openPanels.constraints}
			<div class="panel-content">
				<!-- Barycenter constraint -->
				<div class="constraint-section">
					<label class="checkbox-label">
						<input
							type="checkbox"
							bind:checked={useBarycenterConstraint}
							on:change={updateOptimizationSettings}
						/>
						Fix Barycenter
					</label>
					{#if useBarycenterConstraint}
						<div class="constraint-params">
							<div class="input-group">
								<label>
									X: <input
										type="number"
										bind:value={barycenterTarget[0]}
										on:change={updateOptimizationSettings}
										min="0"
										max="1000"
									/>
								</label>
								<label>
									Y: <input
										type="number"
										bind:value={barycenterTarget[1]}
										on:change={updateOptimizationSettings}
										min="0"
										max="1000"
									/>
								</label>
							</div>
						</div>
					{/if}
				</div>

				<!-- Length constraint -->
				<div class="constraint-section">
					<label class="checkbox-label">
						<input
							type="checkbox"
							bind:checked={useLengthConstraint}
							on:change={updateOptimizationSettings}
						/>
						Fix Curve Length
					</label>
					{#if useLengthConstraint}
						<div class="constraint-params">
							<div class="radio-group">
								<label class="radio-label">
									<input
										type="radio"
										bind:group={useLengthPercentage}
										value={true}
										on:change={updateOptimizationSettings}
									/>
									Percentage of initial length
								</label>
								{#if useLengthPercentage}
									<div class="slider-with-value">
										<input
											type="range"
											bind:value={lengthPercentage}
											on:change={updateOptimizationSettings}
											min="10"
											max="200"
											step="1"
										/>
										<span
											>{lengthPercentage}% of {formatNumber($initialTotalLength)} = {formatNumber(
												($initialTotalLength * lengthPercentage) / 100
											)}</span
										>
									</div>
								{/if}

								<label class="radio-label">
									<input
										type="radio"
										bind:group={useLengthPercentage}
										value={false}
										on:change={updateOptimizationSettings}
									/>
									Absolute length
								</label>
								{#if !useLengthPercentage}
									<div class="input-with-label">
										<input
											type="number"
											bind:value={lengthTarget}
											on:change={updateOptimizationSettings}
											min="0"
											step="1"
										/>
										<span>Current initial length: {formatNumber($initialTotalLength)}</span>
									</div>
								{/if}
							</div>
						</div>
					{/if}
				</div>

				<!-- Constraint implementation method -->
				<div class="constraint-section">
					<label class="checkbox-label">
						<input
							type="checkbox"
							bind:checked={useFullConstraintProjection}
							on:change={updateOptimizationSettings}
						/>
						Use full constraint projection (paper method)
					</label>
				</div>
			</div>
		{/if}
	</div>

	<!-- NUMERICAL PARAMETERS SECTION -->
	<div class="control-panel">
		<div class="panel-header" on:click={() => togglePanel('numerical')}>
			<h4>Numerical Parameters</h4>
			<span class="toggle-icon">{openPanels.numerical ? '▼' : '►'}</span>
		</div>

		{#if openPanels.numerical}
			<div class="panel-content">
				<!-- Line Search Parameters -->
				<div class="control-group">
					<h5>Line Search Parameters</h5>

					<div class="slider-group">
						<label>
							Initial Step Size:
							<div class="slider-with-value">
								<input
									type="range"
									bind:value={$config.tauInitial}
									on:input={updateConfig}
									min="0.01"
									max="1.0"
									step="0.01"
								/>
								<span>{$config.tauInitial.toFixed(2)}</span>
							</div>
						</label>

						<label>
							Armijo Constant (a):
							<div class="slider-with-value">
								<input
									type="range"
									bind:value={$config.aConst}
									on:input={updateConfig}
									min="0.001"
									max="0.5"
									step="0.001"
								/>
								<span>{$config.aConst.toFixed(3)}</span>
							</div>
						</label>

						<label>
							Backtracking Factor (b):
							<div class="slider-with-value">
								<input
									type="range"
									bind:value={$config.bConst}
									on:input={updateConfig}
									min="0.1"
									max="0.9"
									step="0.01"
								/>
								<span>{$config.bConst.toFixed(2)}</span>
							</div>
						</label>

						<label>
							Max Line Search Iterations:
							<div class="slider-with-value">
								<input
									type="range"
									bind:value={$config.maxLineSearch}
									on:input={updateConfig}
									min="10"
									max="100"
									step="1"
								/>
								<span>{$config.maxLineSearch}</span>
							</div>
						</label>
					</div>
				</div>

				<!-- Epsilon Parameters -->
				<div class="control-group">
					<h5>Epsilon Parameters</h5>

					<div class="slider-group">
						<label>
							Epsilon Stability:
							<div class="slider-with-value">
								<input
									type="range"
									bind:value={epsilonStabilityExponent}
									on:input={updateConfig}
									min="-10"
									max="-5"
									step="1"
								/>
								<span>{formatExponent(epsilonStabilityExponent)}</span>
							</div>
						</label>

						<label>
							Epsilon Kernel:
							<div class="slider-with-value">
								<input
									type="range"
									bind:value={epsilonKernelExponent}
									on:input={updateConfig}
									min="-10"
									max="-5"
									step="1"
								/>
								<span>{formatExponent(epsilonKernelExponent)}</span>
							</div>
						</label>

						<label>
							Finite Difference H:
							<div class="slider-with-value">
								<input
									type="range"
									bind:value={finiteDiffHExponent}
									on:input={updateConfig}
									min="-6"
									max="-3"
									step="1"
								/>
								<span>{formatExponent(finiteDiffHExponent)}</span>
							</div>
						</label>

						<label>
							Constraint Tolerance:
							<div class="slider-with-value">
								<input
									type="range"
									bind:value={constraintToleranceExponent}
									on:input={updateConfig}
									min="-10"
									max="-1"
									step="1"
								/>
								<span>{formatExponent(constraintToleranceExponent)}</span>
							</div>
						</label>
					</div>
				</div>
			</div>
		{/if}
	</div>

	<!-- ADVANCED SECTION -->
	<div class="control-panel">
		<div class="panel-header" on:click={() => togglePanel('advanced')}>
			<h4>Advanced Settings</h4>
			<span class="toggle-icon">{openPanels.advanced ? '▼' : '►'}</span>
		</div>

		{#if openPanels.advanced}
			<div class="panel-content">
				<!-- Constraint Stabilization Parameters -->
				<div class="control-group">
					<h5>Constraint Stabilization</h5>

					<div class="slider-group">
						<label>
							Barycenter Stabilization:
							<div class="slider-with-value">
								<input
									type="range"
									bind:value={barycenterStabilization}
									on:input={updateOptimizationSettings}
									min="0.001"
									max="0.1"
									step="0.001"
								/>
								<span>{barycenterStabilization.toFixed(3)}</span>
							</div>
						</label>

						<label>
							Length Stabilization:
							<div class="slider-with-value">
								<input
									type="range"
									bind:value={lengthStabilization}
									on:input={updateOptimizationSettings}
									min="0.001"
									max="0.1"
									step="0.001"
								/>
								<span>{lengthStabilization.toFixed(3)}</span>
							</div>
						</label>
					</div>
				</div>

				<!-- Other advanced settings could go here -->
				<div class="control-group">
					<h5>Other Settings</h5>

					<label class="checkbox-label">
						<input
							type="checkbox"
							bind:checked={$config.applyPerturbation}
							on:change={updateConfig}
						/>
						Apply random perturbation when stuck
					</label>

					<label class="checkbox-label">
						<input
							type="checkbox"
							bind:checked={$config.useSubverticesInEnergy}
							on:change={updateConfig}
						/>
						Include subvertices in energy calculation
					</label>
				</div>
			</div>
		{/if}
	</div>
</div>

<style>
	.controls-container {
		background-color: #f8f9fa;
		border-radius: 8px;
		padding: 15px;
		margin-bottom: 15px;
		width: 350px;
		box-shadow: 0 2px 5px rgba(0, 0, 0, 0.1);
		font-family:
			-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
	}

	.controls-header {
		display: flex;
		justify-content: space-between;
		align-items: center;
		margin-bottom: 15px;
		border-bottom: 1px solid #dee2e6;
		padding-bottom: 8px;
	}

	.user-info {
		display: flex;
		flex-direction: column;
		align-items: flex-end;
		font-size: 0.8rem;
		color: #6c757d;
	}

	h3 {
		margin: 0;
		color: #343a40;
		font-weight: 600;
		font-size: 1.2rem;
	}

	h4 {
		margin: 0;
		font-size: 1.1rem;
		font-weight: 500;
		color: #495057;
	}

	h5 {
		margin: 8px 0;
		font-size: 0.95rem;
		font-weight: 500;
		color: #495057;
	}

	.control-panel {
		margin-bottom: 10px;
		border: 1px solid #dee2e6;
		border-radius: 5px;
		overflow: hidden;
	}

	.panel-header {
		display: flex;
		justify-content: space-between;
		align-items: center;
		background-color: #e9ecef;
		padding: 8px 12px;
		cursor: pointer;
		user-select: none;
		transition: background-color 0.2s;
	}

	.panel-header:hover {
		background-color: #dee2e6;
	}

	.toggle-icon {
		font-size: 0.8rem;
		font-weight: bold;
	}

	.panel-content {
		padding: 12px;
		background-color: #ffffff;
	}

	.control-group {
		margin-bottom: 15px;
		padding-bottom: 10px;
		border-bottom: 1px dashed #dee2e6;
	}

	.control-group:last-child {
		margin-bottom: 0;
		padding-bottom: 0;
		border-bottom: none;
	}

	.constraint-section {
		margin-bottom: 10px;
	}

	.constraint-params {
		margin-left: 25px;
		margin-top: 5px;
		padding: 8px;
		background-color: #f1f3f5;
		border-radius: 4px;
	}

	.input-group {
		display: flex;
		gap: 10px;
		align-items: center;
	}

	input[type='number'] {
		width: 70px;
		padding: 4px;
		border: 1px solid #ced4da;
		border-radius: 4px;
	}

	input[type='range'] {
		width: 100%;
		margin: 8px 0;
	}

	.slider-group {
		display: flex;
		flex-direction: column;
		gap: 12px;
		margin-top: 8px;
	}

	.slider-with-value {
		display: flex;
		align-items: center;
		gap: 10px;
	}

	.slider-with-value span {
		min-width: 85px;
		font-size: 0.85rem;
		color: #495057;
	}

	.input-with-label {
		display: flex;
		align-items: center;
		gap: 10px;
		margin-top: 5px;
		margin-left: 20px;
	}

	.input-with-label span {
		font-size: 0.8rem;
		color: #495057;
	}

	.checkbox-label,
	.radio-label {
		display: flex;
		align-items: center;
		gap: 6px;
		margin-bottom: 6px;
		font-size: 0.9rem;
	}

	.radio-group {
		display: flex;
		flex-direction: column;
		gap: 8px;
	}

	label {
		display: block;
		margin-bottom: 4px;
		font-size: 0.9rem;
		color: #495057;
	}

	input[type='checkbox'],
	input[type='radio'] {
		margin: 0;
	}
</style>
