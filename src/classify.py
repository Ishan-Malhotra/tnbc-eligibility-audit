import asyncio
import json
import os
import aiohttp

# Configuration parameters
CONCURRENT_REQUESTS_LIMIT = 5  # Safe traffic limit for rate ceilings
INPUT_DATA_PATH = 'data/trials.json'
OUTPUT_DATA_PATH = 'data/classified_trials.json'

# Fetch the Anthropic API key from your Mac's environment variables
# (This defaults to None for now, which is perfect while we are just testing the loop)
ANTHROPIC_API_KEY = os.environ.get("ANTHROPIC_API_KEY")

# systemic prompt instructing Claude how to evaluate medical language based on your philosophy
CLASSIFIER_SYSTEM_PROMPT = """
You are an expert molecular oncologist and clinical trial auditor. Your task is to analyze clinical trial eligibility criteria for structural bias against the Basal-Like Immune Suppressed (BLIS) molecular subtype of Triple-Negative Breast Cancer.

Review the trial title, interventions, and eligibility text provided. Determine if the criteria implicitly exclude or under-represent BLIS patients (who present with immune-desert tumors, low TILs, low PD-L1, and rapid recurrence/progression characteristics).

You must respond STRICTLY with a valid JSON object matching this schema:
{
  "nctId": "string",
  "is_immunotherapy": true/false,
  "blia_compatible": true/false,
  "blis_compatible": true/false,
  "exclusionary_mechanisms": ["list", "of", "phrases"],
  "short_biological_rationale": "string structural explanation"
}
"""

async def classify_single_trial(session, semaphore, trial):
    """Handles the async API call to Anthropic for an individual trial record."""
    async with semaphore:
        # Construct the user prompt payload
        user_message = f"Title: {trial['title']}\nInterventions: {json.dumps(trial['interventions'])}\nCriteria:\n{trial['eligibilityCriteria']}"
        
        # This is where your Anthropic API call will go once credits activate
        # For now, we stub the structural loop framework
        try:
            # When credits drop, replace this block with the real aiohttp anthropic payload request
            await asyncio.sleep(0.1) # Simulate quick network passage
            
            # Simulated output matching the target schema structure
            return {
                "nctId": trial['nctId'],
                "is_immunotherapy": any(inv.get('name','').lower() in ['pembrolizumab', 'keytruda', 'pd-1', 'pd-l1'] for inv in trial['interventions']),
                "blia_compatible": True,
                "blis_compatible": False if "PD-L1" in trial['eligibilityCriteria'] else True,
                "exclusionary_mechanisms": ["Biomarker threshold requirements"] if "PD-L1" in trial['eligibilityCriteria'] else [],
                "short_biological_rationale": "Stub validation placeholder"
            }
        except Exception as e:
            print(f"❌ Error processing trial {trial['nctId']}: {str(e)}")
            return None

async def run_batch_classification():
    # Load your fresh 1,452 trial records
    with open(INPUT_DATA_PATH, 'r') as f:
        trials = json.load(f)
        
    print(f"Loaded {len(trials)} trials from local storage. Preparing classification matrix...")
    
    semaphore = asyncio.Semaphore(CONCURRENT_REQUESTS_LIMIT)
    
    # Optional: For testing right now, just slice the first 5 records so it runs instantly
    test_slice = trials[:5] 
    
    async with aiohttp.ClientSession() as session:
        tasks = [classify_single_trial(session, semaphore, trial) for trial in test_slice]
        results = await asyncio.gather(*tasks)
        
    # Filter out empty errors
    final_classified_dataset = [r for r in results if r is not None]
    
    with open(OUTPUT_DATA_PATH, 'w') as f:
        json.dump(final_classified_dataset, f, indent=2)
        
    print(f"✅ Saved classification mapping matrix to {OUTPUT_DATA_PATH}")

if __name__ == '__main__':
    asyncio.run(run_batch_classification())