import requests
import json
import time

def fetch_all_tnbc_trials():
    url = 'https://clinicaltrials.gov/api/v2/studies'
    
    # Pass the raw string. requests will safely handle the URL encoding under the hood!
    params = {
        'query.cond': 'Triple Negative Breast Cancer',
        'countTotal': 'true',
        'pageSize': 100
    }
    
    all_studies = []
    next_page_token = None
    
    print("🚀 Connecting to ClinicalTrials.gov API v2...")
    print(f"🔍 Searching for condition: {params['query.cond']}")
    
    while True:
        if next_page_token:
            params['pageToken'] = next_page_token
            
        response = requests.get(url, params=params)
        
        if response.status_code != 200:
            print(f"❌ Error: Server responded with status {response.status_code}")
            break
            
        data = response.json()
        studies = data.get('studies', [])
        
        if not studies and len(all_studies) == 0:
            print("⚠️ The API returned an empty list. Double-checking response structure...")
            break
            
        all_studies.extend(studies)
        print(f"📥 Pulled {len(all_studies)} / {data.get('totalCount', 'unknown')} trials...")
        
        # Check if another page exists via the token parameter
        next_page_token = data.get('nextPageToken')
        if not next_page_token:
            print("🏁 Reached the final page of available records.")
            break
            
        time.sleep(0.5)
        
    if all_studies:
        cleaned_trials = []
        for study in all_studies:
            protocol = study.get('protocolSection', {})
            ident = protocol.get('identificationModule', {})
            eligibility = protocol.get('eligibilityModule', {})
            arms_module = protocol.get('armsInterventionsModule', {})
            status_mod = protocol.get('statusModule', {})
            
            cleaned_trials.append({
                'nctId': ident.get('nctId'),
                'title': ident.get('officialTitle') or ident.get('briefTitle'),
                'eligibilityCriteria': eligibility.get('eligibilityCriteria', ''),
                'phase': status_mod.get('phases', []),  # API v2 tracks this as a list named 'phases'
                'interventions': arms_module.get('interventions', [])
            })
            
        with open('data/trials.json', 'w') as f:
            json.dump(cleaned_trials, f, indent=2)
            
        print(f"✅ Success! Saved {len(cleaned_trials)} clean trials to data/trials.json")
    else:
        print("❌ No data was collected.")

if __name__ == '__main__':
    fetch_all_tnbc_trials()