from selenium import webdriver

import pytesseract
import cv2
import enzyme
import urllib.request

class get_enzymes:
    def __init__(ge):
        ge.driver = webdriver.Chrome()

    # Compatability with only genscript archive: https://web.archive.org/web/20090502113246/http://www.genscript.com/product_001/enzyme/op/all_ez/start/A/list.html
    # Also references wikipedia for sequences and cut sites

    wikipedia = 'https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites'

    def open_website(ge, website = 'https://www.neb.com/products/restriction-endonucleases'):
    	ge.driver.get(website)


    	for category in driver.find_elements_by_class_name('category'):
    		enzymes = category.find_elements_by_class_name('gae-racking-prod-listing')

    		if len(enzymes) == 0:
    			continue

    		for enzyme in enzymes:
    			enzyme.click()
    			name = driver.find_element_by_class_name('product-detail__title').text
    			
    			sizes = driver.find_elements_by_class_name('product-info__size')

    			prices = driver.find_elements_by_class_name('product-info__listprice')

    			if (sizes[2] >= sizes[1]):
    				size = sizes[2]
    				price = prices[2]
    			else:
    				size = sizes[1]
    				price = prices[1]

    			ge.driver.get(wikipedia)
    			
				navigation = driver.find_elements_by_class_name('toc plainlinks hlist')    			

				

    			driver.get(website)


