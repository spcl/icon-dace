require 'mkmf'
require 'rbconfig'
if RUBY_VERSION[0,3] == '1.8'
  puts "-I#{Config::expand(CONFIG['archdir'])}"
elsif RbConfig::CONFIG['MAJOR'].to_i >= 2
  puts "-I#{RbConfig::CONFIG['rubyhdrdir']} -I#{RbConfig::CONFIG['rubyarchhdrdir']}"
else
  puts "-I#{RbConfig::CONFIG['rubyhdrdir']} -I#{RbConfig::CONFIG['rubyhdrdir']}/#{RbConfig::CONFIG['arch']}"
end
