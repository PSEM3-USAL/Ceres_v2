from django import forms
from django.core.exceptions import ValidationError
from django.utils.translation import ugettext_lazy as _

#class UserInputForm(forms.Form):
    #budget = forms.FloatField(help_text="Enter a budget in USD")
    
    ## Validation (positive number)
    #def clean_budget(self):
        #data = self.cleaned_data['budget'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        
        ##Check if budget is positive
        #if data < 0:
            #raise ValidationError(_('Invalid budget. It must be a positive value')) # wraps this text in one of Django's translation functions ugettext_lazy() (imported as _()), which is good practice if you want to translate your site later.
        
        #return data
        
from django.forms import ModelForm
from cereslibrary.models import UserInputData, AD_UserInputData, Economic_UserInputData

class AD_UserInputForm(ModelForm):
    
    def clean_biogas_product(self):
        biogas_product = self.cleaned_data['biogas_product'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        AD_decision = self.cleaned_data['AD_decision']
        if AD_decision :
            pass
        else:
            if biogas_product != 'None': # discount_rate <= 0:# | discount_rate == None:
                msg = forms.ValidationError("If 'AD' not selected this field must be null")
                self.add_error('biogas_product', msg)
        return biogas_product
    
    class Meta:
        model = AD_UserInputData
        #fields = '__all__'
        fields = ['AD_decision', 'biogas_product']
        
class Economic_UserInputForm(ModelForm):
    
    def clean_discount_rate(self):
        customized_discount_rate = self.cleaned_data.get('customized_discount_rate')
        discount_rate = self.cleaned_data['discount_rate']
        # conditional field 

        if customized_discount_rate == 'cus':
            if discount_rate == None or discount_rate <= 0 : # discount_rate <= 0:# | discount_rate == None:
                msg = forms.ValidationError("This field is required. Only positive values")
                self.add_error('discount_rate', msg)
        elif customized_discount_rate == 'def':
            if discount_rate != None: # discount_rate <= 0:# | discount_rate == None:
                msg = forms.ValidationError("If 'Default' selected this field must be null")
                self.add_error('discount_rate', msg)
            discount_rate = 7
        return discount_rate
    
    def clean_facility_size(self):
        REC_size_clean = self.cleaned_data['REC'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        
        #Check if budget is positive
        if REC_size_clean < 0:
            raise ValidationError(_('Invalid facility size. It must be a positive value')) # wraps this text in one of Django's translation functions ugettext_lazy() (imported as _()), which is good practice if you want to translate your site later.
        
        return REC
    
    def clean_facility_size(self):
        RIN_size_clean = self.cleaned_data['RIN'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        
        #Check if budget is positive
        if RIN_size_clean < 0:
            raise ValidationError(_('Invalid facility size. It must be a positive value')) # wraps this text in one of Django's translation functions ugettext_lazy() (imported as _()), which is good practice if you want to translate your site later.
        
        return RIN
    
    def clean_facility_size(self):
        Capital_Cost_incentive_size_clean = self.cleaned_data['Capital_Cost_incentive'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        
        #Check if budget is positive
        if Capital_Cost_incentive_size_clean < 0:
            raise ValidationError(_('Invalid facility size. It must be a positive value')) # wraps this text in one of Django's translation functions ugettext_lazy() (imported as _()), which is good practice if you want to translate your site later.
        
        return Capital_Cost_incentive
    
    class Meta:
        model = Economic_UserInputData
        #fields = '__all__'
        fields = ['customized_discount_rate', 'discount_rate', 'REC', 'RIN', 'Capital_Cost_incentive']

class UserInputForm(ModelForm):
    """
    # Budget validation (positive number)
    def clean_user_budget(self):
        data = self.cleaned_data['user_budget'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        
        #Check if budget is positive
        if data < 0:
            raise ValidationError(_('Invalid budget. It must be a positive value')) # wraps this text in one of Django's translation functions ugettext_lazy() (imported as _()), which is good practice if you want to translate your site later.
        
        return data
    """
        # Composition validation (positive number)
    #def clean_manure_composition(self):
        #data = self.cleaned_data['manure_composition'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        
        ##Check if budget is positive
        #if data < 0:
            #raise ValidationError(_('Invalid composition. It must be a positive value')) # wraps this text in one of Django's translation functions ugettext_lazy() (imported as _()), which is good practice if you want to translate your site later.
        
        #return data
    
    def clean_number_DairyCow(self):
        number_DairyCow_clean = self.cleaned_data['number_DairyCow'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        
        #Check if budget is positive
        if number_DairyCow_clean < 0:
            raise ValidationError(_('Invalid facility size. It must be a positive value')) # wraps this text in one of Django's translation functions ugettext_lazy() (imported as _()), which is good practice if you want to translate your site later.
        
        return number_DairyCow_clean
    
    def clean_number_DairyHeifer(self):
        number_DairyHeifer_clean = self.cleaned_data['number_DairyHeifer'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        
        #Check if budget is positive
        if number_DairyHeifer_clean < 0:
            raise ValidationError(_('Invalid facility size. It must be a positive value')) # wraps this text in one of Django's translation functions ugettext_lazy() (imported as _()), which is good practice if you want to translate your site later.
        
        return number_DairyHeifer_clean
    
    def clean_number_DairyCalf(self):
        number_DairyCalf_clean = self.cleaned_data['number_DairyCalf'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        
        #Check if budget is positive
        if number_DairyCalf_clean < 0:
            raise ValidationError(_('Invalid facility size. It must be a positive value')) # wraps this text in one of Django's translation functions ugettext_lazy() (imported as _()), which is good practice if you want to translate your site later.
        
        return number_DairyCalf_clean
    
    def clean_number_BeefCow(self):
        number_BeefCow_clean = self.cleaned_data['number_BeefCow'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        
        #Check if budget is positive
        if number_BeefCow_clean < 0:
            raise ValidationError(_('Invalid facility size. It must be a positive value')) # wraps this text in one of Django's translation functions ugettext_lazy() (imported as _()), which is good practice if you want to translate your site later.
        
        return number_BeefCow_clean
    
    def clean_number_BeefCalf(self):
        number_BeefCalf_clean = self.cleaned_data['number_BeefCalf'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        
        #Check if budget is positive
        if number_BeefCalf_clean < 0:
            raise ValidationError(_('Invalid facility size. It must be a positive value')) # wraps this text in one of Django's translation functions ugettext_lazy() (imported as _()), which is good practice if you want to translate your site later.
        
        return number_BeefCalf_clean
    
    def clean_longitude(self):
        longitude_clean = self.cleaned_data['longitude'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        
        #Check if budget is positive
        if longitude_clean > 0:
            raise ValidationError(_('For facilities located in the U.S. the longitude must be negative')) # wraps this text in one of Django's translation functions ugettext_lazy() (imported as _()), which is good practice if you want to translate your site later.
        
        return longitude_clean
    
    def clean_latitude(self):
        latitude_clean = self.cleaned_data['latitude'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        
        #Check if budget is positive
        if latitude_clean < 0:
            raise ValidationError(_('For facilities located in the U.S. the longitude must be positive')) # wraps this text in one of Django's translation functions ugettext_lazy() (imported as _()), which is good practice if you want to translate your site later.
        
        return latitude_clean
    
    #def clean_discharge_mode(self):
        #discharge_mode_clean = self.cleaned_data['discharge_mode'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
               
        #return discharge_mode_clean
    
    #def fields_required(self, fields):
        #"""Used for conditionally marking fields as required."""
        #for field in fields:
            #if not self.cleaned_data.get(field, ''):
                #msg = forms.ValidationError("This field is required.")
                #self.add_error(field, msg)
            
    #def clean_discount_rate(self):
        #custom_discount_rate = self.cleaned_data['custom_discount_rate']

        #if custom_discount_rate:
            #self.fields_required(['discount_rate'])
            #discount_rate_clean = self.cleaned_data['discount_rate']
            ##Check if budget is positive
            ##if discount_rate <= 0:
                ##raise ValidationError(_('Invalid discount rate value. It must be a positive value'))
        #else:
            ##discount_rate_clean = self.cleaned_data['discount_rate'] = 7
            ##self.cleaned_data['discount_rate'] = 7
            #discount_rate_clean = 7

        #return self.cleaned_data, discount_rate_clean
        
        #discount_rate = self.cleaned_data

        return discount_rate
    
    #def clean_AD_decision(self):
        #AD_decision = self.cleaned_data['AD_decision'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
               
        #return AD_decision
    
    #def clean_biogas_product(self):
        #biogas_product = self.cleaned_data['biogas_product'] # This step gets us the data "cleaned" and sanitized of potentially unsafe input using the default validators, and converted into the correct standard type for the data 
        #AD_decision = self.cleaned_data['AD_decision']
        #if AD_decision :
            #pass
        #else:
            #if biogas_product != 'None': # discount_rate <= 0:# | discount_rate == None:
                #msg = forms.ValidationError("If 'AD' not selected this field must be null")
                #self.add_error('biogas_product', msg)
        #return biogas_product
    
    #DEF, CUS = 'def', 'cus'
    
    #SEL = (
        #(DEF, 'Default'),
        #(CUS, 'Custom'),
    #)
    
    #custom_discount_rate = forms.ModelChoiceField(max_length=3,
                            #choices=SEL, widget=forms.RadioSelect)
    
    #discount_rate = forms.FloatField(blank=True, null=True, help_text="(%). Only required if 'Custom discount rate' is selected.")
    
    
    
    class Meta:
        model = UserInputData
        #fields = '__all__'
        #fields = ['facility_size', 'longitude', 'latitude', 'manure_composition','AD_decision', 'discharge_mode', 'customized_discount_rate', 'discount_rate', 'biogas_product']
        fields = ['number_DairyCow', 'number_DairyHeifer', 'number_DairyCalf', 'number_BeefCow', 'number_BeefCalf', 'longitude', 'latitude', 'manure_composition',]
        
    #def __init__(self, *args, **kwargs):
        ##super(UserInputForm, self).__init__(*args, **kwargs)
        #super().__init__(*args, **kwargs)
        ##self.fields['discount_rate'] = custom_discount_rate
        #self.fields['biogas_product'].queryset = custom_discount_rate

        #instance = getattr(self, 'custom_discount_rate', None)
        #if instance == 'def':
            #self.fields['discount_rate'].widget.attrs['readonly'] = True
        
    #def __init__(self, data=None, *args, **kwargs):
        #super(UserInputForm, self).__init__(data, *args, **kwargs)

        ## If 'Custom' is chosen, set send_date as required
        #if data and data.get('custom_discount_rate', None) == self.CUS:
            #self.fields['discount_rate'].required = True
        
    #def __init__(self, *args, **kwargs):
        #super().__init__(*args, **kwargs)
        #self.fields['discount_rate_custom'].queryset = discount_rate_custom.objects.none()
