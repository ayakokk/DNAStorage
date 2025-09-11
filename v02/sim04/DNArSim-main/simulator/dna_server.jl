#---------------------------------------------------------------------------#
#	DNA Channel Server - Persistent Julia Process for C++ Integration     #
#	Modified from original DNArSim simulator for pipe communication	    #
#---------------------------------------------------------------------------#

using DelimitedFiles
include("functions.jl")
include("channel.jl")

# Load probability data for k=4 (investigating default probability issue)
k = 4
# Change to correct directory before loading
cd(dirname(@__FILE__))
include("loadProb.jl")

# Signal that server is ready for communication
println("DNA_SERVER_READY")
flush(stdout)

# Main server loop - process DNA sequences from stdin
while true
    try
        line = strip(readline(stdin))
        # Exit command
        if line == "EXIT" || line == "QUIT"
            break
        end
        
        # Skip empty lines
        if isempty(line)
            continue
        end
        
        # Process DNA sequence
       # Handle string, array, or even nested array inputs
        line_str = line
        while isa(line_str, AbstractArray) && !isempty(line_str)
            line_str = first(line_str)
        end
        # SubStringを完全にクリーンなStringに変換する
        clean_str = String(line_str)
    
        # Now, check if the result is a valid DNA string
        # if文の条件式もclean_strを使うように変更する
        if Base.occursin(r"^[ACGT]+$", clean_str::String)
            # =======================
            # Valid DNA sequence
        # =======================

            # Create sequence array as expected by channel function
            # Redirect Julia's stdout to stderr temporarily to suppress debug messages
            old_stdout = stdout
            redirect_stdout(stderr)
            
            try
                # Run single simulation (nbrSim=1, k=4)
                simSeq = channel(1, k, clean_str)
                
                # Restore stdout
                redirect_stdout(old_stdout)
                
                # Extract and return first simulated sequence
                if !isempty(simSeq) && !isempty(simSeq[1])
                    result_seq = join(simSeq[1])
                    println("RESULT:", result_seq)
                else
                    println("ERROR:NO_SIMULATION_OUTPUT")
                end
            catch e
                # Restore stdout in case of error
                redirect_stdout(old_stdout)
                println("ERROR:JULIA_EXCEPTION:", e)
            end
        else
            println("ERROR:INVALID_DNA_SEQUENCE")
        end
        
        flush(stdout)
        
    catch e
        println("ERROR:JULIA_EXCEPTION:", e)
        flush(stdout)
    end
end

# Clean exit - handle potential broken pipe gracefully
try
    println("DNA_SERVER_EXIT")
    flush(stdout)
catch
    # Ignore all errors during exit (broken pipe is expected)
    # C++ may have already closed the connection
end