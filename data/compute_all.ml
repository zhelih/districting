(* from ExtLib *)
let split str sep =
  let p = String.index str sep in
  let len = 1 in
  let slen = String.length str in
  String.sub str 0 p, String.sub str (p + len) (slen - p - len)

let nsplit str sep =
  if str = "" then []
  else
    let rec loop acc pos =
      if pos > String.length str then
        List.rev acc
      else
        let i = try String.index_from str pos sep with Not_found -> String.length str in
        loop (String.sub str pos (i - pos) :: acc) (i + 1)
    in
    loop [] 0

open Printf
let () =
  (* compute graph from <state>_features.csv *)
  if Array.length Sys.argv < 2 then begin
    printf "Usage: %s <state>\n" Sys.argv.(0);
    exit 0
  end;
  let state = Sys.argv.(1) in

  let f_in = try open_in (state ^ "_features.csv") with _ -> (printf "Input not found\n"; exit 1) in
  (* format: geoid, pop, nbs *)
  (* skip first row *)
  ignore (input_line f_in);

  (* can be done in one read, but I am toooo lazy *)
  let h_pop = Hashtbl.create 1 in
  let h_ind = Hashtbl.create 1 in
  let rec loop i =
    try
      let line = input_line f_in in
(*       printf "1: Processing line %S\n" line; *)
      let (geoid,rest) = split line ',' in
      let (pop,rest) = split rest ',' in
      let pop = try int_of_string pop with _ -> 0 in
      Hashtbl.add h_pop geoid pop;
      Hashtbl.add h_ind geoid i;
      loop (i+1)
    with End_of_file -> ()
  in
  loop 1;
  seek_in f_in 0;

(*   printf "first loop done\n"; *)

  let n = Hashtbl.length h_ind in (* nr_nodes *)
  let adj = Array.make_matrix (n+1) (n+1) false in
  (* skip first row *)
  ignore (input_line f_in);
  let rec loop i =
    try
      let line = input_line f_in in
(*       printf "2: Processing line %S\n" line; *)
      let (geoid,rest) = split line ',' in (* too lazy to merge with code above *)
      let (_,rest) = split rest ',' in
      let rest = if rest.[0] = '"' then String.sub rest 1 (String.length rest - 2) else rest in (* remove quotes *)
      let nbs = nsplit rest ',' in

      let v_from = try Hashtbl.find h_ind geoid with exn -> (printf "missing %s (v_from)\n" geoid; raise exn) in
      List.iter (fun geoid_to ->
        let v_to = try Hashtbl.find h_ind geoid_to with exn -> (printf "missing %s from %s (v_to)\n" geoid_to geoid; raise exn) in
        adj.(v_from).(v_to) <- true;
        adj.(v_to).(v_from) <- true
      ) nbs;
      loop (i+1)
    with End_of_file -> ()
  in
  loop 1;
  close_in f_in;

  let f = open_out (state ^ ".dimacs") in

  fprintf f "p edge %d\n" n;
  for i = 1 to n do
    for j = i+1 to n do
      if adj.(i).(j) then
        fprintf f "e %d %d\n" (i-1) (j-1)
    done
  done;
  close_out f;

  printf "adj done\n";

  let f = open_out (state ^ ".hash") in
  fprintf f "#   Used to map graph vertices to ids\n";
  Hashtbl.iter (fun k v ->
    fprintf f "%d %s\n" (v-1) k
  ) h_ind;
  close_out f;
  let f = open_out (state ^ ".population") in
  (* compute total population *)
  let total_pop = ref 0 in
  Hashtbl.iter (fun _k v -> total_pop := !total_pop + v) h_pop;
  fprintf f "#   Used to get population for a graph vertex\n";
  fprintf f "total pop = %d\n" !total_pop;
  Hashtbl.iter (fun k v ->
    let ind = Hashtbl.find h_ind k in
    fprintf f "%d %d\n" (ind-1) v
  ) h_pop;
  close_out f
