require 'CdoCdiBuild'
class CdiBuild < CdoCdiBuild
  def generateTasks
    super
    @help[:checkV] = 'run "cdo -V" to x-check the build'
    @builders.each {|builder|
      task "#{builder.name}_checkV".to_sym => ["#{builder.name}_make".to_sym] do |t|
        CommandExecution.execute("./app/cdi -V",builder)
      end
    }
  end
end

build = CdiBuild.new(ENV['HOME']+'/.rake_cdi.json')
build.generateTasks
