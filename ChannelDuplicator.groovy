/*
 * meta-sweeper - for performing parametric sweeps of simulated
 * metagenomic sequencing experiments.
 * Copyright (C) 2016 "Matthew Z DeMaere"
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
import groovyx.gpars.dataflow.DataflowQueue

class ChannelDuplicator {
    DataflowQueue orig

    ChannelDuplicator(DataflowQueue orig) {
        this.orig = orig
    }

    DataflowQueue onCopy() {
	    def copied, keep
        (copied, keep) = this.orig.into(2)
        this.orig = keep
        return copied
    }

    static ChannelDuplicator createFrom(Object[] o) {
        return new ChannelDuplicator(Channel.from(o))
    }

    static ChannelDuplicator createFrom(DataflowQueue q) {
        return new ChannelDuplicator(q)
    }
}
